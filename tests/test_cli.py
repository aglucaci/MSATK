from pathlib import Path

from msatk.cli import main

DATA = Path(__file__).parent / "data"


def test_profile_cli(tmp_path):
    code = main(
        ["profile", str(DATA / "dna.fasta"), "--out", str(tmp_path), "--tables-only", "--force"]
    )
    assert code == 0
    assert (tmp_path / "tables" / "alignment_summary.csv").exists()
    assert (tmp_path / "summary.json").exists()


def test_qc_cli(tmp_path):
    code = main(
        ["qc", str(DATA / "dna.fasta"), "--out", str(tmp_path), "--max-gap-site", "0.1", "--force"]
    )
    assert code == 0
    assert (tmp_path / "qc_warnings.txt").exists()


def test_embed_cli(tmp_path):
    code = main(["embed", str(DATA / "dna.fasta"), "--out", str(tmp_path), "--force"])
    assert code == 0
    assert (tmp_path / "sequence_embeddings.csv").exists()


def test_demo_cli(tmp_path):
    code = main(["demo", "--out", str(tmp_path), "--force"])
    assert code == 0
    assert (tmp_path / "report.html").exists()
