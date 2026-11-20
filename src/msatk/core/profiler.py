"""High-level alignment profiler shared by the CLI, Python API and reports."""

from __future__ import annotations

import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from importlib import resources
from pathlib import Path
from typing import Any

from msatk.core.detect import normalize_requested_type
from msatk.core.interpretation import interpret_summary, recommended_next_steps
from msatk.core.qc import QCFinding, QCThresholds, flagged_sequences, flagged_sites, qc_findings
from msatk.io import write_csv, write_json, write_matrix_csv
from msatk.io.readers import parse_partition_file, read_sequence_alignment
from msatk.stats import StatisticsBundle, compute_statistics, gap_block_rows
from msatk.stats.registry import metric_dictionary
from msatk.version import __version__

OUTPUT_SCHEMA_VERSION = "2.0"

PRESETS: dict[str, dict[str, Any]] = {
    # Fast: no all-vs-all comparison, no extra figures.
    "basic": {"pairwise": "none", "alignment_figure": False, "plots": "png"},
    # Default: size-aware pairwise, standard figures, alignment overview.
    "standard": {"pairwise": "auto", "alignment_figure": True, "plots": "png"},
    # Everything, including exact pairwise, PCA embedding and SVG copies of figures.
    "comprehensive": {
        "pairwise": "exact",
        "alignment_figure": True,
        "plots": "png,svg",
        "embeddings": "pca",
        "encodings": ["categorical", "onehot"],
    },
}


class MSATK:
    """Unified MSATK profiling engine for CLI and Python usage.

    Statistics are computed lazily, once, and cached. ``pairwise`` controls the
    cost of all-vs-all comparison: ``auto`` (exact up to ``max_exact_sequences``,
    sampled above), ``exact``, ``sampled`` or ``none``.
    """

    def __init__(
        self,
        path: str | Path,
        sequence_type: str = "auto",
        fmt: str = "auto",
        validation_mode: str = "permissive",
        pairwise: str = "auto",
        max_exact_sequences: int = 1000,
        sample_pairs: int = 20000,
        seed: int = 0,
        genetic_code: int = 1,
        partitions: str | Path | None = None,
    ) -> None:
        self.path = Path(path)
        self.genetic_code = int(genetic_code)
        self.partitions_file = Path(partitions) if partitions else None
        self.validation_mode = validation_mode
        self.alignment, self.annotations = read_sequence_alignment(
            self.path, fmt=fmt, validation_mode=validation_mode
        )
        self.partitions = (
            parse_partition_file(self.partitions_file, self.alignment.length)
            if self.partitions_file
            else []
        )
        self._special: tuple[Any, ...] | None = None
        self.sequence_type = normalize_requested_type(sequence_type, self.alignment)
        self.pairwise_method = pairwise
        self.max_exact_sequences = max_exact_sequences
        self.sample_pairs = sample_pairs
        self.seed = seed
        self._bundle: StatisticsBundle | None = None

    # ------------------------------------------------------------------ statistics
    def stats(self) -> StatisticsBundle:
        if self._bundle is None:
            self._bundle = compute_statistics(
                self.alignment,
                self.sequence_type,
                pairwise=self.pairwise_method,
                max_exact_sequences=self.max_exact_sequences,
                sample_pairs=self.sample_pairs,
                seed=self.seed,
                genetic_code=self.genetic_code,
            )
        return self._bundle

    def summary(self) -> dict[str, Any]:
        return dict(self.stats().summary)

    def sites(self) -> object:
        """Per-site table (DataFrame when pandas is installed, else list of dicts)."""

        return self._frame(self.stats().sites)

    def sequences(self) -> object:
        """Per-sequence table."""

        return self._frame(self.stats().sequences)

    def pairwise(self) -> object:
        """Long-format pairwise table (requires exact pairwise mode)."""

        result = self.stats().pairwise
        if result is None or not result.is_exact:
            self._force_exact_pairwise()
            result = self.stats().pairwise
        return self._frame(result.long_table() if result else [])

    def gap_blocks(self) -> object:
        return self._frame(gap_block_rows(self.stats().encoded))

    def qc(self, thresholds: QCThresholds | None = None) -> list[dict[str, Any]]:
        return [f.as_dict() for f in self.qc_findings(thresholds)]

    def qc_findings(self, thresholds: QCThresholds | None = None) -> list[QCFinding]:
        bundle = self.stats()
        return qc_findings(
            self.alignment,
            self.sequence_type,
            bundle.sequences,
            bundle.sites,
            summary=bundle.summary,
            thresholds=thresholds,
            genetic_code=self.genetic_code,
        )

    # Backward-compatible names -------------------------------------------------
    def per_sequence_stats(self) -> object:
        return self.sequences()

    def per_site_stats(self) -> object:
        return self.sites()

    def composition_summary(self) -> object:
        return self._frame(self.stats().composition)

    def pairwise_identity_matrix(self) -> list[list[float | None]]:
        self._force_exact_pairwise()
        result = self.stats().pairwise
        return result.matrix_rows("identity") if result else [[1.0]]

    def distance_matrix(self) -> list[list[float | None]]:
        self._force_exact_pairwise()
        result = self.stats().pairwise
        return result.matrix_rows("distance") if result else [[0.0]]

    def embeddings(self, method: str = "pca", representation: str = "onehot") -> object:
        from msatk.ml.embeddings import sequence_embeddings

        return self._frame(
            sequence_embeddings(
                self.alignment,
                method=method,
                representation=representation,
                molecule_type=self.sequence_type,
                random_state=self.seed,
            )
        )

    def qc_warnings(self, max_gap_seq: float = 0.3, max_gap_site: float = 0.5) -> list[str]:
        th = QCThresholds(max_gap_seq=max_gap_seq, max_gap_site=max_gap_site)
        return [f.as_text() for f in self.qc_findings(th)]

    def interpretation(self, max_gap_seq: float = 0.3, max_gap_site: float = 0.5) -> str:
        return interpret_summary(self.summary(), self.qc_warnings(max_gap_seq, max_gap_site))

    def recommended_next_steps(
        self, max_gap_seq: float = 0.3, max_gap_site: float = 0.5
    ) -> list[str]:
        return recommended_next_steps(self.summary(), self.qc_warnings(max_gap_seq, max_gap_site))

    # ------------------------------------------------------------------ figures
    def plot_alignment(self, output: str | Path = "alignment.png", **kwargs: Any) -> Path:
        from msatk.viz.alignment import plot_alignment

        return plot_alignment(self.alignment, output, self.sequence_type, **kwargs)

    def plot_logo(self, output: str | Path = "logo.png", **kwargs: Any) -> Path:
        from msatk.viz.logos import plot_logo

        return plot_logo(self.alignment, output, self.sequence_type, **kwargs)

    def plot_entropy(self, path: str | Path = "entropy_by_site.png") -> bool:
        return self._plot_site_series(path, "entropy", "Entropy", "Entropy by site")

    def plot_gap_profile(self, path: str | Path = "gap_fraction_by_site.png") -> bool:
        return self._plot_site_series(path, "gap_fraction", "Gap fraction", "Gap fraction by site")

    def plot_conservation(self, path: str | Path = "conservation_by_site.png") -> bool:
        return self._plot_site_series(
            path, "conservation_score", "Conservation score", "Conservation by site"
        )

    def plot_composition(self, path: str | Path = "composition_summary.png") -> bool:
        from msatk.plots.basic import plot_single_bar

        rows = self.stats().composition
        return plot_single_bar(
            path,
            [str(row["character"]) for row in rows],
            [float(row["fraction"]) for row in rows],
            "Character",
            "Fraction",
            "Character composition",
        )

    # ------------------------------------------------------------------ ML
    def to_tensor(self, encoding: str = "onehot", **kwargs: Any) -> Any:
        """Return an :class:`~msatk.encoding.tensors.EncodedTensor` (``.array`` is NumPy)."""

        from msatk.encoding.tensors import encode

        return encode(self.alignment, self.sequence_type, encoding, **kwargs)

    def to_image(self, output: str | Path = "alignment.png", **kwargs: Any) -> dict[str, Any]:
        from msatk.encoding.export import export_image

        return export_image(self.alignment, self.sequence_type, output, **kwargs)

    def export_encoding(self, out_dir: str | Path, **kwargs: Any) -> dict[str, Any]:
        from msatk.encoding.export import export_tensor

        return export_tensor(self.alignment, self.sequence_type, out_dir, **kwargs)

    # ------------------------------------------------------------------ reports
    def report(
        self, output: str | Path = "report.html", self_contained: bool = True, **kwargs: Any
    ) -> Path:
        """Write an interactive, offline HTML report (single file by default)."""

        from msatk.report.html import render_interactive_report

        summary = self._provenance_summary()
        findings = self.qc_findings()
        target = Path(output)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            render_interactive_report(
                summary,
                self.stats(),
                findings,
                interpretation_text=interpret_summary(summary, [f.as_text() for f in findings]),
                next_steps=recommended_next_steps(summary, [f.as_text() for f in findings]),
                base_dir=target.parent,
                self_contained=self_contained,
                **kwargs,
            ),
            encoding="utf-8",
        )
        return target

    def write_report(
        self, path: str | Path, fmt: str | None = None, figures: list[str] | None = None
    ) -> None:
        target = Path(path)
        fmt = fmt or target.suffix.lstrip(".") or "html"
        if fmt in {"md", "markdown"}:
            from msatk.report.markdown import render_markdown_report

            warnings = self.qc_warnings()
            summary = self.summary()
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(
                render_markdown_report(
                    summary,
                    warnings,
                    figures,
                    interpret_summary(summary, warnings),
                    recommended_next_steps(summary, warnings),
                ),
                encoding="utf-8",
            )
        else:
            self.report(target, figures=figures)

    # ------------------------------------------------------------------ full run
    def write_outputs(
        self,
        outdir: str | Path | None = None,
        tables: bool = True,
        plots: bool = True,
        report: bool = True,
        embeddings_method: str | None = None,
        max_gap_seq: float = 0.3,
        max_gap_site: float = 0.5,
        force: bool = False,
        command: str | None = None,
        plot_formats: str | list[str] = "png",
        alignment_figure: bool = True,
        self_contained: bool = False,
        encodings: list[str] | None = None,
        preset: str | None = None,
    ) -> dict[str, Any]:
        """Run the complete profile and write all outputs to ``outdir``."""

        from msatk.plots.basic import parse_plot_formats

        if preset:
            cfg = PRESETS[preset]
            if self.pairwise_method == "auto" and cfg["pairwise"] != "auto":
                self.pairwise_method = cfg["pairwise"]
                self._bundle = None
            alignment_figure = alignment_figure and cfg["alignment_figure"]
            if plot_formats == "png":
                plot_formats = cfg["plots"]
            embeddings_method = embeddings_method or cfg.get("embeddings")
            encodings = encodings if encodings is not None else cfg.get("encodings")
        formats = parse_plot_formats(plot_formats)
        if not formats:
            plots = False
        out = self._resolve_outdir(outdir, force=force)
        tables_dir = out / "tables"
        figures_dir = out / "figures"
        tables_dir.mkdir(parents=True, exist_ok=True)
        self._write_report_assets(out / "assets")

        bundle = self.stats()
        thresholds = QCThresholds(max_gap_seq=max_gap_seq, max_gap_site=max_gap_site)
        findings = self.qc_findings(thresholds)
        warnings = [f.as_text() for f in findings]
        summary = self._provenance_summary()
        interpretation = interpret_summary(summary, warnings)
        next_steps = recommended_next_steps(summary, warnings)
        summary["interpretation"] = interpretation
        summary["recommended_next_steps"] = next_steps
        embeddings: list[dict[str, Any]] = []
        if embeddings_method and embeddings_method != "none":
            embeddings = self._rows(self.embeddings(method=embeddings_method))
        pw = bundle.pairwise
        figures: list[str] = []
        written: list[str] = []

        if tables:
            write_csv(tables_dir / "alignment_summary.csv", [summary])
            write_csv(tables_dir / "per_sequence_stats.csv", bundle.sequences)
            write_csv(tables_dir / "per_site_stats.csv", bundle.sites)
            write_csv(
                tables_dir / "entropy_by_site.csv",
                [{"site_index": r["site_index"], "entropy": r["entropy"]} for r in bundle.sites],
            )
            write_csv(tables_dir / "composition_summary.csv", bundle.composition)
            write_csv(tables_dir / "gap_blocks.csv", gap_block_rows(bundle.encoded))
            if embeddings:
                write_csv(tables_dir / "sequence_embeddings.csv", embeddings)
            if pw is not None and pw.is_exact:
                write_matrix_csv(
                    tables_dir / "pairwise_identity_matrix.csv",
                    self.alignment.ids,
                    pw.matrix_rows("identity"),
                )
                write_matrix_csv(
                    tables_dir / "distance_matrix.csv",
                    self.alignment.ids,
                    pw.matrix_rows("distance"),
                )
                write_csv(tables_dir / "pairwise_comparisons.csv", pw.long_table())
            write_csv(
                tables_dir / "qc_flags.csv",
                self._qc_flag_rows(bundle.sequences, bundle.sites, max_gap_seq, max_gap_site),
            )
            write_csv(tables_dir / "qc_findings.csv", [f.as_dict() for f in findings])
            write_json(
                out / "qc_findings.json", [f.as_dict() | {"items": f.items} for f in findings]
            )
            write_json(out / "summary.json", summary)
            write_json(out / "metric_dictionary.json", metric_dictionary())
            (out / "qc_warnings.txt").write_text(
                interpretation + "\n\n" + "\n".join(warnings) + "\n", encoding="utf-8"
            )
            self._write_mode_tables(tables_dir)

        if plots:
            from msatk.plots.basic import generate_standard_plots

            figures = generate_standard_plots(
                figures_dir,
                bundle.sequences,
                bundle.sites,
                bundle.composition,
                pw.matrix_rows("identity")
                if pw is not None and pw.identity is not None and len(pw.ids) <= 500
                else None,
                embeddings,
                formats=formats,
                pair_identities=None if pw is None else pw.pair_identities.tolist(),
            )
            if figures and alignment_figure:
                figures += self._alignment_figures(figures_dir, formats)

        if encodings:
            ml_dir = out / "ml"
            for encoding in encodings:
                info = self.export_encoding(ml_dir, encoding=encoding)
                written.append(info["tensor"])

        if report:
            from msatk.report.html import render_interactive_report
            from msatk.report.markdown import render_markdown_report

            (out / "report.html").write_text(
                render_interactive_report(
                    summary,
                    bundle,
                    findings,
                    interpretation_text=interpretation,
                    next_steps=next_steps,
                    figures=figures,
                    base_dir=out,
                    self_contained=self_contained,
                    embeddings=embeddings,
                ),
                encoding="utf-8",
            )
            (out / "report.md").write_text(
                render_markdown_report(summary, warnings, figures, interpretation, next_steps),
                encoding="utf-8",
            )

        config = {
            "tool": "MSATK",
            "version": __version__,
            "input": str(self.path),
            "input_sha256": summary["input_sha256"],
            "sequence_type": self.sequence_type,
            "format": self.alignment.fmt,
            "validation_mode": self.validation_mode,
            "pairwise": self.pairwise_method,
            "max_exact_sequences": self.max_exact_sequences,
            "sample_pairs": self.sample_pairs,
            "seed": self.seed,
            "embeddings": embeddings_method or "none",
            "plot_formats": ",".join(formats) if plots else "none",
            "encodings": ",".join(encodings or []) or "none",
            "preset": preset or "",
            "max_gap_seq": max_gap_seq,
            "max_gap_site": max_gap_site,
            "output_schema_version": OUTPUT_SCHEMA_VERSION,
            "generated_at": summary["generated_at"],
            "python_version": summary["python_version"],
            "platform": summary["platform"],
            "command": command or "",
        }
        (out / "parameters.yaml").write_text(
            "\n".join(f"{k}: {_yaml_scalar(v)}" for k, v in config.items()) + "\n",
            encoding="utf-8",
        )
        write_json(
            out / "run_manifest.json",
            {
                "config": config,
                "dependencies": _dependency_versions(),
                "outputs": sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()),
            },
        )
        (out / "msatk.log").write_text(
            f"MSATK {__version__} profile completed for {self.path}\n"
            f"Detected type: {self.sequence_type}\n"
            f"Pairwise: {summary.get('pairwise_method')}\n"
            f"Findings: {len(findings)}\n",
            encoding="utf-8",
        )
        report_path = out / "report.html"
        return {
            "summary": summary,
            "warnings": warnings,
            "findings": [f.as_dict() for f in findings],
            "interpretation": interpretation,
            "next_steps": next_steps,
            "figures": figures,
            "encodings": written,
            "out_dir": str(out),
            "report": str(report_path) if report and report_path.exists() else "",
        }

    # ------------------------------------------------------------------ helpers
    def _provenance_summary(self) -> dict[str, Any]:
        from msatk.encoding.export import file_sha256

        return {
            **self.summary(),
            "msatk_version": __version__,
            "output_schema_version": OUTPUT_SCHEMA_VERSION,
            "input": str(self.path),
            "input_sha256": file_sha256(self.path),
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "coordinate_convention": "1-based, inclusive",
            "alignment_kind": "sequence",
            "alignment_scope": "single"
            if self.alignment.n_sequences == 1
            else "pairwise"
            if self.alignment.n_sequences == 2
            else "multiple",
            "alignment_method": self.annotations.method,
            "annotations_found": self.annotations.present()
            + (["partitions (file)"] if self.partitions else []),
            "specialized": self.specialized()[0],
            "analysis_applicability": self.specialized()[2],
        }

    def specialized(
        self,
    ) -> tuple[dict[str, dict[str, Any]], dict[str, list[dict[str, Any]]], dict[str, str]]:
        """Analyses that apply only to some inputs (pairwise, protein substitution, RNA
        structure, profile states, partitions): ``(summaries, tables, applicability)``."""

        if self._special is None:
            from msatk.stats.specialized import specialized_analyses

            self._special = specialized_analyses(
                self.stats(), self.sequence_type, self.annotations, self.partitions
            )
        return self._special  # type: ignore[return-value]

    def _force_exact_pairwise(self) -> None:
        result = self.stats().pairwise
        if result is None or not result.is_exact:
            if self.alignment.n_sequences < 2:
                return
            self.pairwise_method = "exact"
            self._bundle = None

    def _plot_site_series(self, path: str | Path, key: str, ylabel: str, title: str) -> bool:
        from msatk.plots.basic import plot_single_series

        rows = self.stats().sites
        return plot_single_series(
            path,
            [int(r["site_index"]) for r in rows],
            [r[key] for r in rows],
            "Site",
            ylabel,
            title,
        )

    def _alignment_figures(self, figures_dir: Path, formats: list[str]) -> list[str]:
        from msatk.viz.alignment import plot_alignment
        from msatk.viz.logos import plot_logo

        names: list[str] = []
        end = min(self.alignment.length, 120)
        try:
            for fmt in formats:
                plot_alignment(
                    self.alignment,
                    figures_dir / f"alignment_overview.{fmt}",
                    self.sequence_type,
                    start=1,
                    end=end,
                    max_rows=60,
                    dpi=200,
                    title=f"Alignment sites 1-{end}"
                    + (
                        f" (first 60 of {self.alignment.n_sequences} sequences)"
                        if self.alignment.n_sequences > 60
                        else ""
                    ),
                )
                plot_logo(
                    self.alignment,
                    figures_dir / f"sequence_logo.{fmt}",
                    self.sequence_type,
                    start=1,
                    end=min(self.alignment.length, 80),
                    dpi=200,
                )
            names += [f"alignment_overview.{formats[0]}", f"sequence_logo.{formats[0]}"]
        except Exception as exc:  # pragma: no cover - figure failures must not abort profiling
            (figures_dir / "ALIGNMENT_FIGURE_ERROR.txt").write_text(str(exc), encoding="utf-8")
        return names

    def _write_mode_tables(self, tables_dir: Path) -> None:
        site_rows = self.stats().sites
        for name, rows in self.specialized()[1].items():
            write_csv(tables_dir / f"{name}.csv", rows)
        if self.sequence_type == "codon":
            from msatk.codon.analysis import (
                amino_acid_from_codons,
                codon_position_entropy,
                codon_position_gc,
                codon_usage,
                rscu,
                stop_codon_report,
            )

            write_csv(
                tables_dir / "codon_usage.csv", codon_usage(self.alignment, self.genetic_code)
            )
            write_csv(tables_dir / "rscu.csv", rscu(self.alignment, self.genetic_code))
            write_csv(tables_dir / "codon_position_gc.csv", codon_position_gc(self.alignment))
            write_csv(tables_dir / "codon_position_entropy.csv", codon_position_entropy(site_rows))
            write_csv(
                tables_dir / "stop_codon_report.csv",
                stop_codon_report(self.alignment, self.genetic_code),
            )
            write_csv(
                tables_dir / "translated_amino_acid_summary.csv",
                amino_acid_from_codons(self.alignment, self.genetic_code),
            )
        elif self.sequence_type in {"protein", "translated_cds"}:
            from msatk.protein.analysis import (
                amino_acid_composition,
                hydrophobicity_summary,
                residue_class_summary,
            )

            write_csv(tables_dir / "amino_acid_summary.csv", amino_acid_composition(self.alignment))
            write_csv(
                tables_dir / "residue_class_summary.csv", residue_class_summary(self.alignment)
            )
            write_csv(
                tables_dir / "hydrophobicity_summary.csv", hydrophobicity_summary(self.alignment)
            )

    def _write_report_assets(self, assets_dir: Path) -> None:
        assets_dir.mkdir(parents=True, exist_ok=True)
        logo = resources.files("msatk.assets").joinpath("msatk_logo.png")
        with resources.as_file(logo) as logo_path:
            shutil.copyfile(logo_path, assets_dir / "msatk_logo.png")

    def _qc_flag_rows(
        self,
        seq_rows: list[dict[str, Any]],
        site_rows: list[dict[str, Any]],
        max_gap_seq: float,
        max_gap_site: float,
    ) -> list[dict[str, Any]]:
        return [
            {**row, "flag_type": "sequence"} for row in flagged_sequences(seq_rows, max_gap_seq)
        ] + [{**row, "flag_type": "site"} for row in flagged_sites(site_rows, max_gap_site)]

    def _resolve_outdir(self, outdir: str | Path | None, force: bool = False) -> Path:
        out = self.path.with_name(f"{self.path.stem}_msatk") if outdir is None else Path(outdir)
        if force or not out.exists():
            return out
        index = 2
        candidate = out.with_name(f"{out.name}_{index}")
        while candidate.exists():
            index += 1
            candidate = out.with_name(f"{out.name}_{index}")
        return candidate

    def _frame(self, rows: list[dict[str, Any]]) -> object:
        try:
            import pandas as pd

            return pd.DataFrame(rows)
        except Exception:
            return rows

    def _rows(self, value: object) -> list[dict[str, Any]]:
        if hasattr(value, "to_dict"):
            import math

            rows = value.to_dict(orient="records")  # type: ignore[call-arg]
            return [
                {k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in r.items()}
                for r in rows
            ]
        return value  # type: ignore[return-value]


def _yaml_scalar(value: object) -> str:
    if value is None or value == "":
        return "''"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (int, float)):
        return str(value)
    return json.dumps(str(value))


def _dependency_versions() -> dict[str, str | None]:
    from importlib.metadata import PackageNotFoundError, version

    out: dict[str, str | None] = {}
    for name in ("numpy", "matplotlib", "pandas", "scikit-learn", "umap-learn", "pysam", "Pillow"):
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = None
    return out


AlignmentProfiler = MSATK
