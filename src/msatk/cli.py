"""MSATK command-line interface."""

from __future__ import annotations

import argparse
import shutil
import sys
from importlib import resources
from pathlib import Path
from typing import Any

from msatk import __version__
from msatk.codon import CodonProfiler
from msatk.core import MSATK
from msatk.io import write_csv
from msatk.protein import ProteinProfiler
from msatk.report import render_html_report, render_markdown_report

EPILOG = """examples:
  msatk profile alignment.fasta
  msatk profile alignment.fasta --out results --type codon --plots png,svg
  msatk qc alignment.fasta --max-gap-seq 0.3 --max-gap-site 0.5
  msatk demo
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="msatk",
        description="MSATK: one command to profile, visualize, and report multiple sequence alignments.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"MSATK {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    profile = subparsers.add_parser(
        "profile", help="Generate the full FastQC-like MSATK report.", epilog=EPILOG
    )
    _add_input_output(profile)
    profile.add_argument(
        "--type",
        default="auto",
        choices=["auto", "dna", "rna", "protein", "codon", "translated_cds", "mixed"],
    )
    profile.add_argument("--format", default="auto", help="Input format; default auto.")
    profile.add_argument(
        "--report", default="html", choices=["html", "none"], help="Report format; default html."
    )
    profile.add_argument(
        "--plots",
        default="png",
        help="Comma-separated plot formats. PNG is currently generated when matplotlib is installed.",
    )
    profile.add_argument(
        "--threads", type=int, default=1, help="Reserved for workflow compatibility."
    )
    profile.add_argument(
        "--strict",
        action="store_true",
        help="Fail on duplicate IDs, unequal lengths, or invalid characters.",
    )
    profile.add_argument(
        "--permissive",
        action="store_true",
        help="Warn and continue where scientifically reasonable.",
    )
    profile.add_argument("--tables-only", action="store_true")
    profile.add_argument("--plots-only", action="store_true")
    profile.add_argument("--embedding-method", default="pca", choices=["pca", "umap", "tsne"])
    profile.add_argument("--max-gap-seq", type=float, default=0.3)
    profile.add_argument("--max-gap-site", type=float, default=0.5)
    profile.add_argument(
        "--force", action="store_true", help="Reuse/overwrite the requested output directory."
    )
    profile.set_defaults(func=cmd_profile)

    qc = subparsers.add_parser("qc", help="Run alignment quality-control tables and warnings.")
    _add_input_output(qc)
    qc.add_argument("--type", default="auto")
    qc.add_argument("--max-gap-seq", type=float, default=0.3)
    qc.add_argument("--max-gap-site", type=float, default=0.5)
    qc.add_argument("--force", action="store_true")
    qc.set_defaults(func=cmd_qc)

    codon = subparsers.add_parser("codon", help="Run codon-aware MSATK analysis.")
    _add_input_output(codon)
    codon.add_argument("--force", action="store_true")
    codon.set_defaults(func=cmd_codon)

    protein = subparsers.add_parser("protein", help="Run protein-specific MSATK analysis.")
    _add_input_output(protein)
    protein.add_argument("--force", action="store_true")
    protein.set_defaults(func=cmd_protein)

    embed = subparsers.add_parser("embed", help="Generate sequence embeddings.")
    _add_input_output(embed)
    embed.add_argument("--method", default="pca", choices=["pca", "umap", "tsne"])
    embed.add_argument("--representation", default="onehot", choices=["onehot", "kmer"])
    embed.add_argument("--force", action="store_true")
    embed.set_defaults(func=cmd_embed)

    report = subparsers.add_parser(
        "report", help="Build a report from an existing MSATK output directory."
    )
    report.add_argument("results_dir")
    report.add_argument("--format", default="html", choices=["html", "md", "markdown"])
    report.set_defaults(func=cmd_report)

    demo = subparsers.add_parser("demo", help="Run MSATK on a bundled example alignment.")
    demo.add_argument("--out", default="msatk_demo", help="Demo output directory.")
    demo.add_argument("--force", action="store_true")
    demo.set_defaults(func=cmd_demo)
    return parser


def _add_input_output(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input alignment file.")
    parser.add_argument(
        "--out", default=None, help="Output directory. Defaults to <input_stem>_msatk."
    )


def cmd_profile(args: argparse.Namespace) -> int:
    try:
        validation_mode = "strict" if args.strict else "permissive"
        profiler = MSATK(
            args.input, sequence_type=args.type, fmt=args.format, validation_mode=validation_mode
        )
        tables = not args.plots_only
        plots = not args.tables_only and args.plots.lower() != "none"
        report = args.report != "none" and not args.tables_only and not args.plots_only
        result = profiler.write_outputs(
            outdir=args.out,
            tables=tables,
            plots=plots,
            report=report,
            embeddings_method=args.embedding_method,
            max_gap_seq=args.max_gap_seq,
            max_gap_site=args.max_gap_site,
            force=args.force,
            command=" ".join(["msatk", *sys.argv[1:]]) if sys.argv else None,
        )
    except Exception as exc:
        _print_error(exc)
        return 2
    _print_done(result)
    return 0


def cmd_qc(args: argparse.Namespace) -> int:
    profiler = MSATK(args.input, sequence_type=args.type)
    out = profiler._resolve_outdir(args.out, force=args.force)
    out.mkdir(parents=True, exist_ok=True)
    seq_rows = profiler._rows(profiler.per_sequence_stats())
    site_rows = profiler._rows(profiler.per_site_stats())
    warnings = profiler.qc_warnings(args.max_gap_seq, args.max_gap_site)
    write_csv(
        out / "qc_flags.csv",
        profiler._qc_flag_rows(seq_rows, site_rows, args.max_gap_seq, args.max_gap_site),
    )
    (out / "qc_warnings.txt").write_text(
        profiler.interpretation(args.max_gap_seq, args.max_gap_site)
        + "\n\n"
        + "\n".join(warnings)
        + "\n",
        encoding="utf-8",
    )
    print(f"MSATK QC complete\nResults written to:\n{out}/")
    return 0


def cmd_codon(args: argparse.Namespace) -> int:
    profiler = CodonProfiler(args.input)
    result = profiler.write_outputs(outdir=args.out, force=args.force)
    _print_done(result)
    return 0


def cmd_protein(args: argparse.Namespace) -> int:
    profiler = ProteinProfiler(args.input)
    result = profiler.write_outputs(outdir=args.out, force=args.force)
    _print_done(result)
    return 0


def cmd_embed(args: argparse.Namespace) -> int:
    profiler = MSATK(args.input)
    out = profiler._resolve_outdir(args.out, force=args.force)
    out.mkdir(parents=True, exist_ok=True)
    write_csv(
        out / "sequence_embeddings.csv",
        profiler._rows(profiler.embeddings(args.method, args.representation)),
    )
    print(f"MSATK embeddings complete\nResults written to:\n{out}/sequence_embeddings.csv")
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    import csv
    import json

    results = Path(args.results_dir)
    summary_path = results / "summary.json"
    warnings_path = results / "qc_warnings.txt"
    figures_dir = results / "figures"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    else:
        with (results / "tables" / "alignment_summary.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            summary = next(csv.DictReader(handle))
    warning_lines = (
        warnings_path.read_text(encoding="utf-8").splitlines() if warnings_path.exists() else []
    )
    warnings = [line for line in warning_lines if line.startswith(("WARNING", "OK"))]
    figures = (
        sorted(path.name for path in figures_dir.glob("*.png")) if figures_dir.exists() else []
    )
    interpretation = str(summary.get("interpretation", ""))
    next_steps = summary.get("recommended_next_steps", [])
    if not isinstance(next_steps, list):
        next_steps = [str(next_steps)]
    if args.format in {"md", "markdown"}:
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
    out = Path(args.out)
    if out.exists() and args.force:
        shutil.rmtree(out)
    source = resources.files("msatk.data.example_alignments").joinpath("demo_codon.fasta")
    with resources.as_file(source) as demo_path:
        profiler = MSATK(demo_path)
        result = profiler.write_outputs(outdir=out, force=args.force, command="msatk demo")
    _print_done(result)
    return 0


def _print_done(result: dict[str, Any]) -> None:
    summary = result["summary"]
    warnings = [warning for warning in result["warnings"] if str(warning).startswith("WARNING")]
    print("\nMSATK alignment profile complete\n")
    print(f"Input: {summary.get('input')}")
    print(f"Detected type: {summary['molecule_type']}")
    print(f"Sequences: {int(summary['number_of_sequences']):,}")
    print(f"Alignment length: {int(summary['alignment_length']):,}")
    print(f"Gap fraction: {float(summary['gap_fraction']):.1%}")
    print(f"Variable sites: {int(summary['variable_sites']):,}")
    print(f"Mean pairwise identity: {float(summary['mean_pairwise_identity']):.1%}")
    print("\nSummary interpretation:")
    print(result["interpretation"])
    print("\nWarnings:")
    if warnings:
        for warning in warnings[:8]:
            print(f"- {warning.removeprefix('WARNING: ').strip()}")
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
