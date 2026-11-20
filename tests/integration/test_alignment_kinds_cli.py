"""CLI coverage for inspect, formats and kind-dispatched profiling."""

from __future__ import annotations

import json

import pytest

from msatk.cli import main


def test_formats_command(capsys):
    assert main(["formats"]) == 0
    assert "| MAF (multiple alignment format) |" in capsys.readouterr().out
    assert main(["formats", "--format", "json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert {"kinds", "formats", "levels"} <= set(data)


@pytest.mark.parametrize(
    "name,kind",
    [("blocks.maf", "block"), ("reads.sam", "read"), ("tiny.gfa", "graph"),
     ("tiny.needle", "sequence"), ("rna_structure.sto", "sequence")],
)  # fmt: skip
def test_inspect_command(data_dir, capsys, name, kind):
    assert main(["inspect", str(data_dir / "formats" / name), "--json"]) == 0
    info = json.loads(capsys.readouterr().out)
    assert info["kind"] == kind and info["status"] == "readable"


@pytest.mark.parametrize(
    "name,table",
    [("blocks.maf", "genomes.csv"), ("reads.sam", "references.csv"), ("tiny.gfa", "segments.csv"),
     ("tiny.paf", "queries.csv"), ("tiny.hmm", "match_states.csv")],
)  # fmt: skip
def test_profile_dispatches_by_kind(data_dir, tmp_path, name, table):
    out = tmp_path / "out"
    assert main(["profile", str(data_dir / "formats" / name), "--out", str(out)]) == 0
    assert (out / "tables" / table).exists()
    summary = json.loads((out / "summary.json").read_text())
    assert summary["alignment_kind"] != "sequence"
    assert (out / "report.html").exists() and (out / "report.md").exists()


def test_profile_sequence_with_specialized_tables(data_dir, tmp_path):
    out = tmp_path / "rna"
    args = ["profile", str(data_dir / "formats" / "rna_structure.sto"), "--out", str(out),
            "--plots", "none", "--no-alignment-figure"]  # fmt: skip
    assert main(args) == 0
    assert (out / "tables" / "rna_base_pairs.csv").exists()
    summary = json.loads((out / "summary.json").read_text())
    assert summary["specialized"]["rna_structure"]["consensus_base_pairs"] == 3
    assert 'id="specialized"' in (out / "report.html").read_text()


def test_genetic_code_flag(data_dir, tmp_path, capsys):
    fasta = tmp_path / "mito.fasta"
    fasta.write_text(">a\nATGTGAAAATAA\n>b\nATGTGGAAATAA\n")
    assert main(["qc", str(fasta), "--type", "codon", "--genetic-code", "2"]) == 0
    assert "INTERNAL_STOP_CODON" not in capsys.readouterr().out
