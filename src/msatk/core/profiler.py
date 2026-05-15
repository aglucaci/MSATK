"""High-level alignment profiler."""

from __future__ import annotations

import platform
import shutil
import sys
from datetime import datetime, timezone
from importlib import resources
from pathlib import Path
from typing import Any

from msatk.core.detect import normalize_requested_type
from msatk.core.interpretation import interpret_summary, recommended_next_steps
from msatk.core.qc import flagged_sequences, flagged_sites, qc_warnings
from msatk.core.stats import (
    alignment_summary,
    composition_summary,
    distance_matrix,
    identity_matrix,
    per_sequence_stats,
    per_site_stats,
)
from msatk.io import read_alignment, write_csv, write_json, write_matrix_csv
from msatk.ml import sequence_embeddings
from msatk.plots import generate_standard_plots
from msatk.plots.basic import plot_single_bar, plot_single_series
from msatk.report import render_html_report, render_markdown_report
from msatk.version import __version__


class MSATK:
    """Unified MSATK profiling engine for CLI and Python usage."""

    def __init__(
        self,
        path: str | Path,
        sequence_type: str = "auto",
        fmt: str = "auto",
        validation_mode: str = "permissive",
    ) -> None:
        self.path = Path(path)
        self.validation_mode = validation_mode
        self.alignment = read_alignment(self.path, fmt=fmt, validation_mode=validation_mode)
        self.sequence_type = normalize_requested_type(sequence_type, self.alignment)

    def summary(self) -> dict[str, Any]:
        return alignment_summary(self.alignment, self.sequence_type)

    def per_sequence_stats(self) -> object:
        return self._frame(per_sequence_stats(self.alignment))

    def per_site_stats(self) -> object:
        return self._frame(per_site_stats(self.alignment))

    def composition_summary(self) -> object:
        return self._frame(composition_summary(self.alignment))

    def pairwise_identity_matrix(self) -> list[list[float]]:
        return identity_matrix(self.alignment)

    def distance_matrix(self) -> list[list[float]]:
        return distance_matrix(self.alignment)

    def embeddings(self, method: str = "pca", representation: str = "onehot") -> object:
        return self._frame(
            sequence_embeddings(self.alignment, method=method, representation=representation)
        )

    def qc_warnings(self, max_gap_seq: float = 0.3, max_gap_site: float = 0.5) -> list[str]:
        return qc_warnings(
            self.alignment,
            self.sequence_type,
            self._rows(self.per_sequence_stats()),
            self._rows(self.per_site_stats()),
            max_gap_seq=max_gap_seq,
            max_gap_site=max_gap_site,
        )

    def interpretation(self, max_gap_seq: float = 0.3, max_gap_site: float = 0.5) -> str:
        warnings = self.qc_warnings(max_gap_seq=max_gap_seq, max_gap_site=max_gap_site)
        return interpret_summary(self.summary(), warnings)

    def recommended_next_steps(
        self, max_gap_seq: float = 0.3, max_gap_site: float = 0.5
    ) -> list[str]:
        warnings = self.qc_warnings(max_gap_seq=max_gap_seq, max_gap_site=max_gap_site)
        return recommended_next_steps(self.summary(), warnings)

    def write_report(
        self, path: str | Path, fmt: str | None = None, figures: list[str] | None = None
    ) -> None:
        target = Path(path)
        fmt = fmt or target.suffix.lstrip(".") or "html"
        warnings = self.qc_warnings()
        interpretation = interpret_summary(self.summary(), warnings)
        next_steps = recommended_next_steps(self.summary(), warnings)
        text = (
            render_markdown_report(self.summary(), warnings, figures, interpretation, next_steps)
            if fmt in {"md", "markdown"}
            else render_html_report(self.summary(), warnings, figures, interpretation, next_steps)
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")

    def plot_entropy(self, path: str | Path = "entropy_by_site.png") -> bool:
        rows = self._rows(self.per_site_stats())
        return plot_single_series(
            path,
            [int(row["site_index"]) for row in rows],
            [float(row["entropy"]) for row in rows],
            "Site",
            "Entropy",
            "MSATK Entropy by Site",
        )

    def plot_gap_profile(self, path: str | Path = "gap_fraction_by_site.png") -> bool:
        rows = self._rows(self.per_site_stats())
        return plot_single_series(
            path,
            [int(row["site_index"]) for row in rows],
            [float(row["gap_fraction"]) for row in rows],
            "Site",
            "Gap fraction",
            "MSATK Gap Fraction by Site",
        )

    def plot_conservation(self, path: str | Path = "conservation_by_site.png") -> bool:
        rows = self._rows(self.per_site_stats())
        return plot_single_series(
            path,
            [int(row["site_index"]) for row in rows],
            [float(row["conservation_score"]) for row in rows],
            "Site",
            "Conservation score",
            "MSATK Conservation by Site",
        )

    def plot_composition(self, path: str | Path = "composition_summary.png") -> bool:
        rows = self._rows(self.composition_summary())
        return plot_single_bar(
            path,
            [str(row["character"]) for row in rows],
            [float(row["fraction"]) for row in rows],
            "Character",
            "Fraction",
            "MSATK Composition Summary",
        )

    def write_outputs(
        self,
        outdir: str | Path | None = None,
        tables: bool = True,
        plots: bool = True,
        report: bool = True,
        embeddings_method: str = "pca",
        max_gap_seq: float = 0.3,
        max_gap_site: float = 0.5,
        force: bool = False,
        command: str | None = None,
    ) -> dict[str, Any]:
        out = self._resolve_outdir(outdir, force=force)
        tables_dir = out / "tables"
        figures_dir = out / "figures"
        assets_dir = out / "assets"
        for directory in (tables_dir, figures_dir):
            directory.mkdir(parents=True, exist_ok=True)
        self._write_report_assets(assets_dir)

        base_summary = self.summary()
        seq_rows = self._rows(self.per_sequence_stats())
        site_rows = self._rows(self.per_site_stats())
        composition = self._rows(self.composition_summary())
        identities = self.pairwise_identity_matrix()
        distances = self.distance_matrix()
        embeddings = self._rows(self.embeddings(method=embeddings_method))
        warnings = qc_warnings(
            self.alignment, self.sequence_type, seq_rows, site_rows, max_gap_seq, max_gap_site
        )
        interpretation = interpret_summary(base_summary, warnings)
        next_steps = recommended_next_steps(base_summary, warnings)
        summary = {
            **base_summary,
            "msatk_version": __version__,
            "output_schema_version": "1.0",
            "input": str(self.path),
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "interpretation": interpretation,
            "recommended_next_steps": next_steps,
        }
        figures: list[str] = []

        if tables:
            write_csv(tables_dir / "alignment_summary.csv", [summary])
            write_csv(tables_dir / "per_sequence_stats.csv", seq_rows)
            write_csv(tables_dir / "per_site_stats.csv", site_rows)
            write_csv(
                tables_dir / "entropy_by_site.csv",
                [{"site_index": row["site_index"], "entropy": row["entropy"]} for row in site_rows],
            )
            write_csv(tables_dir / "composition_summary.csv", composition)
            write_csv(tables_dir / "sequence_embeddings.csv", embeddings)
            write_matrix_csv(
                tables_dir / "pairwise_identity_matrix.csv", self.alignment.ids, identities
            )
            write_matrix_csv(tables_dir / "distance_matrix.csv", self.alignment.ids, distances)
            write_csv(
                tables_dir / "qc_flags.csv",
                self._qc_flag_rows(seq_rows, site_rows, max_gap_seq, max_gap_site),
            )
            write_json(out / "summary.json", summary)
            (out / "qc_warnings.txt").write_text(
                interpretation + "\n\n" + "\n".join(warnings) + "\n", encoding="utf-8"
            )
            self._write_mode_tables(tables_dir, site_rows)

        if plots:
            figures = generate_standard_plots(
                figures_dir, seq_rows, site_rows, composition, identities, embeddings
            )

        if report:
            (out / "report.html").write_text(
                render_html_report(summary, warnings, figures, interpretation, next_steps),
                encoding="utf-8",
            )
            (out / "report.md").write_text(
                render_markdown_report(summary, warnings, figures, interpretation, next_steps),
                encoding="utf-8",
            )
        (out / "parameters.yaml").write_text(
            "\n".join(
                [
                    "tool: MSATK",
                    f"version: {__version__}",
                    f"input: {self.path}",
                    f"sequence_type: {self.sequence_type}",
                    f"format: {self.alignment.fmt}",
                    f"validation_mode: {self.validation_mode}",
                    "output_schema_version: '1.0'",
                    f"generated_at: {summary['generated_at']}",
                    f"python_version: {summary['python_version']}",
                    f"platform: {summary['platform']}",
                    f"command: {command or ''}",
                    f"max_gap_seq: {max_gap_seq}",
                    f"max_gap_site: {max_gap_site}",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        (out / "msatk.log").write_text(
            f"MSATK profile completed for {self.path}\nDetected type: {self.sequence_type}\nWarnings: {len(warnings)}\n",
            encoding="utf-8",
        )
        report_path = out / "report.html"
        return {
            "summary": summary,
            "warnings": warnings,
            "interpretation": interpretation,
            "next_steps": next_steps,
            "figures": figures,
            "out_dir": str(out),
            "report": str(report_path) if report and report_path.exists() else "",
        }

    def _write_mode_tables(self, tables_dir: Path, site_rows: list[dict[str, Any]]) -> None:
        if self.sequence_type == "codon":
            from msatk.codon.analysis import (
                amino_acid_from_codons,
                codon_position_entropy,
                codon_position_gc,
                codon_usage,
                rscu,
                stop_codon_report,
            )

            write_csv(tables_dir / "codon_usage.csv", codon_usage(self.alignment))
            write_csv(tables_dir / "rscu.csv", rscu(self.alignment))
            write_csv(tables_dir / "codon_position_gc.csv", codon_position_gc(self.alignment))
            write_csv(tables_dir / "codon_position_entropy.csv", codon_position_entropy(site_rows))
            write_csv(tables_dir / "stop_codon_report.csv", stop_codon_report(self.alignment))
            write_csv(
                tables_dir / "translated_amino_acid_summary.csv",
                amino_acid_from_codons(self.alignment),
            )
        elif self.sequence_type == "protein":
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
            return value.to_dict(orient="records")  # type: ignore[no-any-return,call-arg]
        return value  # type: ignore[return-value]


AlignmentProfiler = MSATK
