"""Build the MSATK manuscript PDF from Markdown using Pandoc."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "paper" / "paper.md"
DEFAULT_OUTPUT = ROOT / "paper" / "paper.pdf"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert paper/paper.md to a PDF with Pandoc and references.bib."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=DEFAULT_INPUT,
        type=Path,
        help="Markdown manuscript path. Defaults to paper/paper.md.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT,
        type=Path,
        help="Output PDF path. Defaults to paper/paper.pdf.",
    )
    parser.add_argument(
        "--pdf-engine",
        default="xelatex",
        help="Pandoc PDF engine to use, such as xelatex, lualatex, pdflatex, or tectonic.",
    )
    parser.add_argument(
        "--csl",
        type=Path,
        help="Optional CSL style file for reference formatting.",
    )
    return parser.parse_args()


def require_executable(name: str, install_hint: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"Missing required executable: {name}\n{install_hint}")


def main() -> int:
    args = parse_args()
    manuscript = args.input.resolve()
    output = args.output.resolve()

    if not manuscript.exists():
        raise SystemExit(f"Markdown manuscript not found: {manuscript}")

    require_executable(
        "pandoc",
        "Install with: conda install -c conda-forge pandoc",
    )
    require_executable(
        args.pdf_engine,
        "Install a LaTeX engine, for example: conda install -c conda-forge tectonic "
        "and rerun with --pdf-engine tectonic.",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "pandoc",
        str(manuscript),
        "--from",
        "markdown",
        "--to",
        "pdf",
        "--citeproc",
        "--pdf-engine",
        args.pdf_engine,
        "--output",
        str(output),
    ]
    if args.csl:
        command.extend(["--csl", str(args.csl.resolve())])

    subprocess.run(command, cwd=manuscript.parent, check=True)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"PDF build failed with exit code {exc.returncode}", file=sys.stderr)
        raise SystemExit(exc.returncode)
