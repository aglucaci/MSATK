"""Agreement and runtime comparison between MSATK and AMAS (Borowiec 2016).

For each alignment, both tools compute a summary; shared metrics are compared
and wall time / peak memory are measured with ``/usr/bin/time -l`` (macOS) or
``-v`` (Linux) over ``--repeats`` runs (median reported). Each run is a fresh
process, so interpreter start-up is included for both tools.

Usage::

    python paper/analysis/compare_tools.py --amas-python /path/to/venv/bin/python

AMAS must be importable by ``--amas-python`` (``pip install amas``).
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE / "data"
RESULTS = HERE / "results"

REAL = [
    ("PF00042 globin (seed)", "PF00042.fasta", "protein"),
    ("PF00069 protein kinase (seed)", "PF00069.fasta", "protein"),
    ("PF00071 Ras (seed)", "PF00071.fasta", "protein"),
    ("PF00076 RRM (seed)", "PF00076.fasta", "protein"),
    ("H3N2 HA CDS (NCBI)", "h3n2_ha.fasta", "dna"),
]
SYNTHETIC = [(1000, 10000), (5000, 10000)]


def timed(cmd: list[str]) -> tuple[float, float, str]:
    """Run ``cmd`` under /usr/bin/time; return (seconds, peak MB, stdout)."""
    flag = "-l" if sys.platform == "darwin" else "-v"
    proc = subprocess.run(["/usr/bin/time", flag, *cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd} failed:\n{proc.stderr}")
    err = proc.stderr
    if sys.platform == "darwin":
        secs = float(re.search(r"([\d.]+) real", err).group(1))
        rss_mb = int(re.search(r"(\d+)\s+maximum resident set size", err).group(1)) / 2**20
    else:
        m, s = re.search(r"Elapsed \(wall clock\) time.*: (?:(\d+):)?([\d.]+)", err).groups()[-2:]
        secs = float(s) + 60 * float(m or 0)
        rss_mb = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", err).group(1)) / 1024
    return secs, rss_mb, proc.stdout


def run_amas(py: str, path: Path, kind: str, workdir: Path) -> tuple[float, float, dict[str, str]]:
    out = workdir / "amas.tsv"
    secs, mb, _ = timed(
        [
            py,
            "-m",
            "amas.AMAS",
            "summary",
            "-f",
            "fasta",
            "-d",
            "aa" if kind == "protein" else "dna",
            "-i",
            str(path),
            "-o",
            str(out),
        ]
    )
    with out.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    return secs, mb, row


def run_msatk(path: Path, kind: str, workdir: Path) -> tuple[float, float, dict[str, object]]:
    out = workdir / "msatk.json"
    msatk = shutil.which("msatk") or "msatk"
    secs, mb, _ = timed(
        [
            msatk,
            "stats",
            str(path),
            "--type",
            kind,
            "--pairwise",
            "none",
            "--output-format",
            "json",
            "--out",
            str(out),
        ]
    )
    data = json.loads(out.read_text())
    return secs, mb, data[0] if isinstance(data, list) else data


def compare(name: str, path: Path, kind: str, py: str, repeats: int) -> dict[str, object]:
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        amas_runs = [run_amas(py, path, kind, work) for _ in range(repeats)]
        msatk_runs = [run_msatk(path, kind, work) for _ in range(repeats)]
    a, m = amas_runs[0][2], msatk_runs[0][2]
    total = int(m["total_cells"])
    # AMAS "undetermined" = gaps + missing-data symbols (?, N/X); MSATK reports them separately.
    msatk_missing = 100 * (int(m["total_gaps"]) + int(m["missing_cells"])) / total
    row = {
        "dataset": name,
        "type": kind,
        "sequences": m["number_of_sequences"],
        "sites": m["alignment_length"],
        "amas_taxa": int(a["No_of_taxa"]),
        "amas_length": int(a["Alignment_length"]),
        "amas_missing_pct": float(a["Missing_percent"]),
        "msatk_missing_pct": round(msatk_missing, 3),
        "amas_variable": int(a["No_variable_sites"]),
        "msatk_variable": int(m["variable_sites"]),
        "amas_parsimony_informative": int(a["Parsimony_informative_sites"]),
        "msatk_parsimony_informative": int(m["parsimony_informative_sites"]),
        "msatk_ambiguous_cells": int(m["ambiguous_cells"]),
        "amas_seconds": statistics.median(r[0] for r in amas_runs),
        "amas_peak_mb": statistics.median(r[1] for r in amas_runs),
        "msatk_seconds": statistics.median(r[0] for r in msatk_runs),
        "msatk_peak_mb": statistics.median(r[1] for r in msatk_runs),
    }
    row["agree"] = (
        row["amas_taxa"] == row["sequences"]
        and row["amas_length"] == row["sites"]
        and abs(row["amas_missing_pct"] - row["msatk_missing_pct"]) < 0.001
        and row["amas_variable"] == row["msatk_variable"]
        and row["amas_parsimony_informative"] == row["msatk_parsimony_informative"]
    )
    return row


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--amas-python", required=True)
    ap.add_argument("--repeats", type=int, default=3)
    args = ap.parse_args()
    sys.path.insert(0, str(HERE.parent.parent / "benchmarks"))
    from run_benchmarks import synthetic_fasta

    rows = []
    for name, file, kind in REAL:
        rows.append(compare(name, DATA / file, kind, args.amas_python, args.repeats))
        print(rows[-1])
    with tempfile.TemporaryDirectory() as tmp:
        for n, length in SYNTHETIC:
            path = Path(tmp) / f"syn_{n}x{length}.fasta"
            synthetic_fasta(path, n, length)
            rows.append(
                compare(
                    f"Synthetic DNA {n:,} x {length:,}", path, "dna", args.amas_python, args.repeats
                )
            )
            print(rows[-1])
    RESULTS.mkdir(exist_ok=True)
    with (RESULTS / "tool_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    amas_version = subprocess.run(
        [args.amas_python, "-m", "pip", "show", "amas"], capture_output=True, text=True
    ).stdout
    env = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version.split()[0],
        "amas": re.search(r"Version: (\S+)", amas_version).group(1),
        "repeats": args.repeats,
    }
    (RESULTS / "tool_comparison.json").write_text(
        json.dumps({"environment": env, "results": rows}, indent=2)
    )


if __name__ == "__main__":
    main()
