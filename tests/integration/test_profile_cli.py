import json

from msatk.cli import main


def test_profile_tiny_dna(data_dir, tmp_path):
    outdir = tmp_path / "results"
    result = main(
        ["profile", str(data_dir / "dna" / "tiny_dna.fasta"), "--out", str(outdir), "--force"]
    )
    assert result == 0
    assert (outdir / "report.html").exists()
    assert (outdir / "summary.json").exists()
    assert (outdir / "tables" / "alignment_summary.csv").exists()
    assert (outdir / "tables" / "per_sequence_stats.csv").exists()
    assert (outdir / "tables" / "per_site_stats.csv").exists()
    summary = json.loads((outdir / "summary.json").read_text())
    assert summary["number_of_sequences"] == 3


def test_profile_valid_codon_alignment(data_dir, tmp_path):
    outdir = tmp_path / "codon_results"
    result = main(
        [
            "profile",
            str(data_dir / "codon" / "valid_codon_alignment.fasta"),
            "--type",
            "codon",
            "--out",
            str(outdir),
            "--force",
        ]
    )
    assert result == 0
    assert (outdir / "tables" / "codon_usage.csv").exists()
    assert (outdir / "tables" / "rscu.csv").exists()
    assert (outdir / "tables" / "stop_codon_report.csv").exists()
