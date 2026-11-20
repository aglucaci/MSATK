"""MSATK command-line interface."""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import shutil
import sys
from importlib import resources
from pathlib import Path
from typing import Any

from msatk import __version__
from msatk.exceptions import MSATKError

MOLECULE_TYPES = ["auto", "dna", "rna", "protein", "codon", "translated_cds", "mixed"]

EPILOG = """examples:
  msatk profile alignment.fasta --out results/
  msatk profile alignment.fasta --preset comprehensive --plots png,svg,pdf --self-contained
  msatk stats alignment.fasta --level site > sites.csv
  msatk qc alignment.fasta --max-gap-seq 0.3 --max-gap-site 0.5
  msatk plot alignment.fasta --start 1 --end 300 --out region.svg
  msatk logo alignment.fasta --start 1 --end 60 --out logo.pdf
  msatk encode alignment.fasta --encoding onehot --window 256 --stride 128 --out ml/
  msatk image alignment.fasta --kind rgb --scale 4 --out alignment.png
  msatk inspect genomes.maf.gz
  msatk profile reads.bam --out reads_qc/
  msatk profile cds.fasta --type codon --genetic-code 2
  msatk profile supermatrix.phy --partitions partitions.txt
  msatk formats
  msatk demo
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="msatk",
        description="MSATK: statistical profiling, visualization and ML representations of "
        "multiple sequence alignments.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"MSATK {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # profile ----------------------------------------------------------------
    profile = sub.add_parser(
        "profile", help="Full profile: statistics, QC, figures and an offline HTML report.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter,
    )  # fmt: skip
    _add_input(profile)
    _add_out(profile)
    profile.add_argument(
        "--preset", choices=["basic", "standard", "comprehensive"], default=None,
        help="Computation profile. basic: no pairwise matrix; comprehensive: exact pairwise, "
        "PCA embedding, SVG figures and ML tensors.",
    )  # fmt: skip
    profile.add_argument("--report", default="html", choices=["html", "none"])
    profile.add_argument(
        "--self-contained", action="store_true",
        help="Embed figures in report.html so it is a single portable file.",
    )  # fmt: skip
    profile.add_argument(
        "--plots", default="png",
        help="Comma-separated figure formats: png,svg,pdf,tiff or none (default png).",
    )  # fmt: skip
    profile.add_argument(
        "--no-alignment-figure", action="store_true",
        help="Skip the residue-level alignment overview and logo figures.",
    )  # fmt: skip
    _add_pairwise(profile)
    profile.add_argument(
        "--embeddings", default=None, choices=["none", "pca", "umap", "tsne"],
        help="Optional sequence embedding (off by default).",
    )  # fmt: skip
    profile.add_argument(
        "--embedding-method", default=None, choices=["pca", "umap", "tsne"], help=argparse.SUPPRESS
    )
    profile.add_argument(
        "--encodings", default=None,
        help="Comma-separated ML encodings to export under ml/ (categorical,onehot,features,physicochemical).",
    )  # fmt: skip
    profile.add_argument(
        "--threads", type=int, default=1, help="Reserved for workflow compatibility."
    )
    profile.add_argument(
        "--strict",
        action="store_true",
        help="Fail on duplicate IDs, unequal lengths or invalid characters.",
    )
    profile.add_argument("--permissive", action="store_true", help="Warn and continue (default).")
    profile.add_argument("--tables-only", action="store_true")
    profile.add_argument("--plots-only", action="store_true")
    profile.add_argument("--max-gap-seq", type=float, default=0.3)
    profile.add_argument("--max-gap-site", type=float, default=0.5)
    profile.add_argument(
        "--force", action="store_true", help="Reuse/overwrite the output directory."
    )
    profile.add_argument(
        "--reference", default=None, help="Reference FASTA for CRAM input (when required)."
    )
    profile.set_defaults(func=cmd_profile)

    # stats ------------------------------------------------------------------
    stats = sub.add_parser(
        "stats", help="Print statistics at one resolution (CSV/TSV/JSON to stdout or file)."
    )
    _add_input(stats)
    stats.add_argument(
        "--level",
        default="summary",
        choices=["summary", "site", "sequence", "pairwise", "composition", "gaps"],
    )
    stats.add_argument("--output-format", default="csv", choices=["csv", "tsv", "json"])
    stats.add_argument("--out", default=None, help="Output file (default stdout).")
    _add_pairwise(stats)
    stats.set_defaults(func=cmd_stats)

    # qc ---------------------------------------------------------------------
    qc = sub.add_parser("qc", help="Run QC rules and write findings and flags.")
    _add_input(qc)
    _add_out(qc)
    qc.add_argument("--max-gap-seq", type=float, default=0.3)
    qc.add_argument("--max-gap-site", type=float, default=0.5)
    _add_pairwise(qc)
    qc.add_argument("--force", action="store_true")
    qc.set_defaults(func=cmd_qc)

    # plot -------------------------------------------------------------------
    plot = sub.add_parser(
        "plot", help="Render a residue-level alignment figure (PNG/SVG/PDF/TIFF)."
    )
    _add_input(plot)
    plot.add_argument("--out", required=True, help="Output figure path; format from suffix.")
    plot.add_argument("--start", type=int, default=None)
    plot.add_argument("--end", type=int, default=None)
    plot.add_argument("--max-rows", type=int, default=200)
    plot.add_argument("--order", default="input", choices=["input", "consensus_identity", "gaps"])
    plot.add_argument(
        "--theme",
        default="auto",
        help="nucleotide, nucleotide_cvd, protein, grayscale or a JSON/YAML colour file.",
    )
    plot.add_argument(
        "--tracks",
        default="conservation,gaps",
        help="Comma list of conservation,gaps,entropy,occupancy or none.",
    )
    plot.add_argument("--highlight-differences", action="store_true")
    plot.add_argument("--no-letters", action="store_true")
    plot.add_argument("--no-consensus", action="store_true")
    plot.add_argument("--dpi", type=int, default=300)
    plot.add_argument("--title", default=None)
    plot.set_defaults(func=cmd_plot)

    logo = sub.add_parser("logo", help="Render a sequence logo for a region.")
    _add_input(logo)
    logo.add_argument("--out", required=True)
    logo.add_argument("--start", type=int, default=None)
    logo.add_argument("--end", type=int, default=None)
    logo.add_argument("--theme", default="auto")
    logo.add_argument(
        "--no-correction", action="store_true", help="Disable small-sample correction."
    )
    logo.add_argument("--dpi", type=int, default=300)
    logo.set_defaults(func=cmd_logo)

    # encode / image -----------------------------------------------------------
    encode = sub.add_parser("encode", help="Export ML-ready tensors (NPZ) with a manifest.")
    _add_input(encode)
    encode.add_argument("--out", required=True, help="Output directory.")
    encode.add_argument(
        "--encoding",
        default="onehot",
        choices=["categorical", "onehot", "physicochemical", "features"],
    )
    encode.add_argument(
        "--ambiguity", default="missing", choices=["missing", "fractional", "separate"]
    )
    encode.add_argument("--start", type=int, default=None)
    encode.add_argument("--end", type=int, default=None)
    encode.add_argument(
        "--window", type=int, default=None, help="Tile width in sites (no resampling)."
    )
    encode.add_argument("--stride", type=int, default=None)
    encode.add_argument("--rows-per-tile", type=int, default=None)
    encode.add_argument(
        "--label",
        action="append",
        default=[],
        help="KEY=VALUE label stored in the manifest (repeatable).",
    )
    encode.set_defaults(func=cmd_encode)

    image = sub.add_parser("image", help="Export an RGB or lossless label image with a manifest.")
    _add_input(image)
    image.add_argument("--out", required=True, help="Image path (.png, .tif or .npy).")
    image.add_argument("--kind", default="rgb", choices=["rgb", "label"])
    image.add_argument("--theme", default="auto")
    image.add_argument("--scale", type=int, default=1, help="Integer pixels per cell (lossless).")
    image.add_argument(
        "--size", default=None, help="HxW nearest-neighbour resize, e.g. 224x224 (lossy)."
    )
    image.add_argument("--start", type=int, default=None)
    image.add_argument("--end", type=int, default=None)
    image.set_defaults(func=cmd_image)

    # legacy specialised commands ------------------------------------------------
    codon = sub.add_parser("codon", help="Run codon-aware MSATK analysis.")
    codon.add_argument("input")
    _add_out(codon)
    codon.add_argument("--force", action="store_true")
    codon.set_defaults(func=cmd_codon)

    protein = sub.add_parser("protein", help="Run protein-specific MSATK analysis.")
    protein.add_argument("input")
    _add_out(protein)
    protein.add_argument("--force", action="store_true")
    protein.set_defaults(func=cmd_protein)

    embed = sub.add_parser("embed", help="Generate 2-D sequence embeddings.")
    _add_input(embed)
    _add_out(embed)
    embed.add_argument("--method", default="pca", choices=["pca", "umap", "tsne"])
    embed.add_argument("--representation", default="onehot", choices=["onehot", "kmer"])
    embed.add_argument("--force", action="store_true")
    embed.set_defaults(func=cmd_embed)

    report = sub.add_parser(
        "report", help="Rebuild report.html from an existing MSATK output directory."
    )
    report.add_argument("results_dir")
    report.add_argument("--format", default="html", choices=["html", "md", "markdown"])
    report.add_argument("--self-contained", action="store_true")
    report.set_defaults(func=cmd_report)

    demo = sub.add_parser("demo", help="Run MSATK on a bundled example alignment.")
    demo.add_argument("--out", default="msatk_demo")
    demo.add_argument("--force", action="store_true")
    demo.set_defaults(func=cmd_demo)

    metrics = sub.add_parser("metrics", help="Print the metric dictionary (Markdown or JSON).")
    metrics.add_argument("--format", default="markdown", choices=["markdown", "json"])
    metrics.set_defaults(func=cmd_metrics)

    inspect = sub.add_parser(
        "inspect",
        help="Report what an input contains, its data model and which analyses are valid.",
    )
    inspect.add_argument("input")
    inspect.add_argument("--format", default="auto")
    inspect.add_argument("--json", action="store_true", help="Print JSON.")
    inspect.set_defaults(func=cmd_inspect)

    formats = sub.add_parser(
        "formats", help="Print the format/alignment-type compatibility matrix."
    )
    formats.add_argument("--format", default="markdown", choices=["markdown", "json"])
    formats.set_defaults(func=cmd_formats)
    return parser


def _add_input(p: argparse.ArgumentParser) -> None:
    p.add_argument("input", help="Input alignment file.")
    p.add_argument("--type", default="auto", choices=MOLECULE_TYPES)
    p.add_argument("--format", default="auto", help="Input format; default auto.")
    p.add_argument(
        "--genetic-code", type=int, default=1,
        help="NCBI translation table for codon/stop analysis (default 1, standard).",
    )  # fmt: skip
    p.add_argument(
        "--partitions", default=None,
        help="Partition file (RAxML/IQ-TREE style or NEXUS charsets) for per-partition stats.",
    )  # fmt: skip


def _add_out(p: argparse.ArgumentParser) -> None:
    p.add_argument("--out", default=None, help="Output directory. Defaults to <input_stem>_msatk.")


def _add_pairwise(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--pairwise", default="auto", choices=["auto", "exact", "sampled", "none"],
        help="Pairwise identity: auto = exact up to --max-exact-sequences, sampled above.",
    )  # fmt: skip
    p.add_argument("--max-exact-sequences", type=int, default=1000)
    p.add_argument("--sample-pairs", type=int, default=20000)
    p.add_argument("--seed", type=int, default=0)


def _profiler(args: argparse.Namespace, cls: Any = None, **extra: Any) -> Any:
    from msatk.core import MSATK

    kwargs: dict[str, Any] = {
        "sequence_type": getattr(args, "type", "auto"),
        "fmt": getattr(args, "format", "auto"),
    }
    for name in (
        "pairwise", "max_exact_sequences", "sample_pairs", "seed", "genetic_code", "partitions"
    ):  # fmt: skip
        if hasattr(args, name):
            kwargs[name] = getattr(args, name)
    kwargs.update(extra)
    return (cls or MSATK)(args.input, **kwargs)


def cmd_profile(args: argparse.Namespace) -> int:
    try:
        from msatk.alignments import detect_input_format, get_format, load_alignment
        from msatk.alignments.outputs import write_kind_outputs

        fmt = args.format if args.format != "auto" else detect_input_format(args.input)
        kind = get_format(fmt).kind if fmt in {f.name for f in _formats()} else "sequence"
        if kind != "sequence":
            obj = load_alignment(args.input, fmt, reference=args.reference)
            res = write_kind_outputs(
                obj, args.input, args.out, force=args.force,
                command=" ".join(["msatk", *sys.argv[1:]]) if sys.argv else None,
            )  # fmt: skip
            print(f"\nMSATK {kind} profile complete\n")
            for key, value in list(res["summary"].items())[:14]:
                print(f"  {key}: {value}")
            print(f"\nOutputs: {res['outdir']}")
            return 0
        profiler = _profiler(args, validation_mode="strict" if args.strict else "permissive")
        result = profiler.write_outputs(
            outdir=args.out,
            tables=not args.plots_only,
            plots=not args.tables_only and args.plots.lower() != "none",
            report=args.report != "none" and not args.tables_only and not args.plots_only,
            embeddings_method=args.embeddings or args.embedding_method,
            max_gap_seq=args.max_gap_seq,
            max_gap_site=args.max_gap_site,
            force=args.force,
            command=" ".join(["msatk", *sys.argv[1:]]) if sys.argv else None,
            plot_formats=args.plots if args.plots.lower() != "none" else "png",
            alignment_figure=not args.no_alignment_figure,
            self_contained=args.self_contained,
            encodings=[e for e in args.encodings.split(",") if e] if args.encodings else None,
            preset=args.preset,
        )
    except (MSATKError, ValueError, OSError) as exc:
        _print_error(exc)
        return 2
    _print_done(result)
    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    try:
        profiler = _profiler(args)
        bundle = profiler.stats()
        if args.level == "summary":
            rows: list[dict[str, Any]] | dict[str, Any] = bundle.summary
        elif args.level == "site":
            rows = bundle.sites
        elif args.level == "sequence":
            rows = bundle.sequences
        elif args.level == "composition":
            rows = bundle.composition
        elif args.level == "gaps":
            from msatk.stats import gap_block_rows

            rows = gap_block_rows(bundle.encoded)
        else:
            if args.pairwise in {"sampled", "none"}:
                raise ValueError("--level pairwise needs --pairwise exact or auto.")
            profiler._force_exact_pairwise()
            pw = profiler.stats().pairwise
            rows = pw.long_table() if pw else []
    except (MSATKError, ValueError, OSError) as exc:
        _print_error(exc)
        return 2
    with (
        open(args.out, "w", encoding="utf-8", newline="")
        if args.out
        else contextlib.nullcontext(sys.stdout)
    ) as handle:
        _write_rows(handle, rows, args.output_format)
    return 0


def _write_rows(handle: Any, rows: list[dict[str, Any]] | dict[str, Any], fmt: str) -> None:
    if fmt == "json":
        json.dump(rows, handle, indent=2, default=str)
        handle.write("\n")
        return
    delimiter = "\t" if fmt == "tsv" else ","
    if isinstance(rows, dict):
        writer = csv.writer(handle, delimiter=delimiter)
        writer.writerow(["statistic", "value"])
        for key, value in rows.items():
            writer.writerow([key, "" if value is None else value])
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    dict_writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
    dict_writer.writeheader()
    dict_writer.writerows(rows)


def cmd_qc(args: argparse.Namespace) -> int:
    from msatk.core.qc import QCThresholds
    from msatk.io import write_csv, write_json

    profiler = _profiler(args)
    out = profiler._resolve_outdir(args.out, force=args.force)
    out.mkdir(parents=True, exist_ok=True)
    bundle = profiler.stats()
    findings = profiler.qc_findings(QCThresholds(args.max_gap_seq, args.max_gap_site))
    write_csv(
        out / "qc_flags.csv",
        profiler._qc_flag_rows(bundle.sequences, bundle.sites, args.max_gap_seq, args.max_gap_site),
    )
    write_csv(out / "qc_findings.csv", [f.as_dict() for f in findings])
    write_json(out / "qc_findings.json", [f.as_dict() | {"items": f.items} for f in findings])
    (out / "qc_warnings.txt").write_text(
        profiler.interpretation(args.max_gap_seq, args.max_gap_site)
        + "\n\n"
        + "\n".join(f.as_text() for f in findings)
        + "\n",
        encoding="utf-8",
    )
    for f in findings:
        print(f"[{f.severity.upper():7}] {f.code}: {f.message}")
    print(f"\nMSATK QC complete\nResults written to:\n{out}/")
    return 0


def cmd_plot(args: argparse.Namespace) -> int:
    try:
        tracks = () if args.tracks == "none" else tuple(t for t in args.tracks.split(",") if t)
        path = _profiler(args, pairwise="none").plot_alignment(
            args.out,
            start=args.start,
            end=args.end,
            max_rows=args.max_rows,
            order=args.order,
            theme=args.theme,
            tracks=tracks,
            highlight_differences=args.highlight_differences,
            show_letters=False if args.no_letters else None,
            show_consensus=not args.no_consensus,
            dpi=args.dpi,
            title=args.title,
        )
    except (MSATKError, ValueError, OSError) as exc:
        _print_error(exc)
        return 2
    print(f"MSATK alignment figure written: {path}")
    return 0


def cmd_logo(args: argparse.Namespace) -> int:
    try:
        path = _profiler(args, pairwise="none").plot_logo(
            args.out,
            start=args.start,
            end=args.end,
            theme=args.theme,
            dpi=args.dpi,
            small_sample_correction=not args.no_correction,
        )
    except (MSATKError, ValueError, OSError) as exc:
        _print_error(exc)
        return 2
    print(f"MSATK sequence logo written: {path}")
    return 0


def cmd_encode(args: argparse.Namespace) -> int:
    labels = dict(item.split("=", 1) for item in args.label if "=" in item)
    try:
        info = _profiler(args, pairwise="none").export_encoding(
            args.out,
            encoding=args.encoding,
            ambiguity=args.ambiguity,
            start=args.start,
            end=args.end,
            window=args.window,
            stride=args.stride,
            rows_per_tile=args.rows_per_tile,
            labels=labels,
        )
    except (MSATKError, ValueError, OSError) as exc:
        _print_error(exc)
        return 2
    print(f"MSATK tensor written: {info['tensor']} shape={tuple(info['shape'])}")
    print(f"Manifest: {info['manifest']}")
    return 0


def cmd_image(args: argparse.Namespace) -> int:
    size = None
    if args.size:
        h, w = args.size.lower().split("x")
        size = (int(h), int(w))
    try:
        info = _profiler(args, pairwise="none").to_image(
            args.out, kind=args.kind, theme=args.theme, scale=args.scale, size=size,
            start=args.start, end=args.end,
        )  # fmt: skip
    except (MSATKError, ValueError, OSError) as exc:
        _print_error(exc)
        return 2
    print(f"MSATK image written: {info['image']} shape={tuple(info['shape'])}")
    if size:
        print("Note: resizing is lossy; pixel coordinates no longer map 1:1 to alignment cells.")
    return 0


def cmd_codon(args: argparse.Namespace) -> int:
    from msatk.codon import CodonProfiler

    result = CodonProfiler(args.input).write_outputs(outdir=args.out, force=args.force)
    _print_done(result)
    return 0


def cmd_protein(args: argparse.Namespace) -> int:
    from msatk.protein import ProteinProfiler

    result = ProteinProfiler(args.input).write_outputs(outdir=args.out, force=args.force)
    _print_done(result)
    return 0


def cmd_embed(args: argparse.Namespace) -> int:
    from msatk.io import write_csv

    try:
        profiler = _profiler(args, pairwise="none")
        out = profiler._resolve_outdir(args.out, force=args.force)
        out.mkdir(parents=True, exist_ok=True)
        rows = profiler._rows(profiler.embeddings(args.method, args.representation))
    except (MSATKError, ValueError, OSError) as exc:
        _print_error(exc)
        return 2
    write_csv(out / "sequence_embeddings.csv", rows)
    print(
        f"MSATK embeddings ({rows[0]['method'] if rows else args.method}) written to:\n{out}/sequence_embeddings.csv"
    )
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    from msatk.report import render_html_report, render_markdown_report

    results = Path(args.results_dir)
    summary_path = results / "summary.json"
    figures_dir = results / "figures"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    else:
        with (results / "tables" / "alignment_summary.csv").open(encoding="utf-8", newline="") as h:
            summary = next(csv.DictReader(h))
    figures = (
        sorted(p.name for p in figures_dir.iterdir() if p.suffix in {".png", ".svg"})
        if figures_dir.exists()
        else []
    )
    interpretation = str(summary.get("interpretation", ""))
    next_steps = summary.get("recommended_next_steps", [])
    if not isinstance(next_steps, list):
        next_steps = [str(next_steps)]
    warnings_path = results / "qc_warnings.txt"
    warnings = [
        line
        for line in (
            warnings_path.read_text(encoding="utf-8").splitlines() if warnings_path.exists() else []
        )
        if line.startswith(("WARNING", "OK", "NOTE", "ERROR"))
    ]
    source = Path(str(summary.get("input", "")))
    if args.format == "html" and source.is_file():
        # Re-render the full interactive report from the original input.
        from msatk.core import MSATK

        profiler = MSATK(source, sequence_type=str(summary.get("molecule_type", "auto")))
        target = profiler.report(
            results / "report.html", self_contained=args.self_contained, figures=figures
        )
    elif args.format in {"md", "markdown"}:
        target = results / "report.md"
        target.write_text(
            render_markdown_report(summary, warnings, figures, interpretation, next_steps),
            encoding="utf-8",
        )
    else:
        target = results / "report.html"
        target.write_text(
            render_html_report(summary, warnings, figures, interpretation, next_steps),
            encoding="utf-8",
        )
    print(f"MSATK report written: {target}")
    return 0


def cmd_demo(args: argparse.Namespace) -> int:
    from msatk.core import MSATK

    out = Path(args.out)
    if out.exists() and args.force:
        shutil.rmtree(out)
    source = resources.files("msatk.data.example_alignments").joinpath("demo_codon.fasta")
    with resources.as_file(source) as demo_path:
        result = MSATK(demo_path).write_outputs(
            outdir=out, force=args.force, command="msatk demo", self_contained=True
        )
    _print_done(result)
    return 0


def cmd_metrics(args: argparse.Namespace) -> int:
    from msatk.stats.registry import metric_dictionary, metrics_markdown

    if args.format == "json":
        print(json.dumps(metric_dictionary(), indent=2))
    else:
        print(metrics_markdown())
    return 0


def _pct(value: Any) -> str:
    return "NA" if value is None else f"{float(value):.1%}"


def _formats() -> Any:
    from msatk.alignments import FORMATS

    return FORMATS


def cmd_inspect(args: argparse.Namespace) -> int:
    from msatk.alignments import inspect_alignment

    try:
        info = inspect_alignment(args.input, args.format)
    except (MSATKError, ValueError, OSError, KeyError) as exc:
        _print_error(exc)
        return 2
    if args.json:
        print(json.dumps(info, indent=2, default=str))
        return 0
    print(f"\n{info['path']}")
    print(f"  format:        {info['format_label']} ({info['format']}, {info['support_level']})")
    if info.get("compressed"):
        print(f"  compressed:    {info['compressed']}")
    print(f"  data model:    {info['data_model']} (kind: {info['kind']})")
    for key in ("molecule_type", "scope", "method", "sequences", "columns"):
        if key in info:
            print(f"  {key + ':':<14} {info[key]}")
    if info.get("annotations_found"):
        print(f"  annotations:   {', '.join(info['annotations_found'])}")
    for key, value in (info.get("contents") or {}).items():
        if value is not None:
            print(f"  {key + ':':<34} {value}")
    print("\n  valid analyses:")
    for item in info["valid_analyses"]:
        print(f"    - {item}")
    if info.get("specialized_analyses"):
        print("\n  specialized analyses for this input:")
        for name, why in info["specialized_analyses"].items():
            print(f"    - {name}: {why}")
    if info["not_computed_by_design"]:
        print("\n  not computed by design:")
        for item in info["not_computed_by_design"]:
            print(f"    - {item}")
    print(f"\n  status: {info['status']}\n")
    return 0


def cmd_formats(args: argparse.Namespace) -> int:
    from msatk.alignments import capability_registry, compatibility_markdown

    if args.format == "json":
        print(json.dumps(capability_registry(), indent=2))
    else:
        print(compatibility_markdown())
    return 0


def _print_done(result: dict[str, Any]) -> None:
    summary = result["summary"]
    warnings = [w for w in result["warnings"] if str(w).startswith(("WARNING", "ERROR"))]
    print("\nMSATK alignment profile complete\n")
    print(f"Input: {summary.get('input')}")
    print(f"Detected type: {summary['molecule_type']}")
    print(f"Sequences: {int(summary['number_of_sequences']):,}")
    print(f"Alignment length: {int(summary['alignment_length']):,}")
    print(f"Gap fraction: {_pct(summary['gap_fraction'])}")
    print(f"Variable sites: {int(summary['variable_sites']):,}")
    print(
        f"Mean pairwise identity: {_pct(summary.get('mean_pairwise_identity'))}"
        f" ({summary.get('pairwise_method', 'none')})"
    )
    print("\nSummary interpretation:")
    print(result["interpretation"])
    print("\nWarnings:")
    if warnings:
        for warning in warnings[:8]:
            print(f"- {warning.split(': ', 1)[-1].strip()}")
        if len(warnings) > 8:
            print(f"- {len(warnings) - 8} additional warning(s); see qc_warnings.txt")
    else:
        print("- No major QC warnings triggered by current thresholds")
    print("\nResults written to:")
    print(f"{result['out_dir']}/")
    if result.get("report"):
        print("Open report:")
        print(result["report"])


def _print_error(exc: Exception) -> None:
    print("\nMSATK could not complete the requested analysis.\n", file=sys.stderr)
    print(str(exc), file=sys.stderr)
    print("\nTry:", file=sys.stderr)
    print("  msatk profile alignment.fasta --type auto", file=sys.stderr)
    print("  msatk demo", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
