from pathlib import Path

from msatk import MSATK, AlignmentProfiler, CodonProfiler, ProteinProfiler, msatk
from msatk.io import read_alignment

DATA = Path(__file__).parent / "data"


def test_read_fasta_and_summary():
    alignment = read_alignment(DATA / "dna.fasta")
    assert alignment.n_sequences == 3
    assert alignment.length == 9
    profiler = MSATK(DATA / "dna.fasta")
    summary = profiler.summary()
    assert summary["tool"] == "MSATK"
    assert summary["number_of_sequences"] == 3
    assert summary["molecule_type"] == "codon"
    assert summary["variable_sites"] >= 1
    assert summary["detected_file_format"] == "fasta"
    assert summary["detected_molecule_type"] == "codon"
    assert summary["sequence_lengths_consistent"] is True
    assert summary["length_divisible_by_three"] is True
    assert summary["stop_codons_exist"] is True
    assert "frameshift_warnings_exist" in summary
    assert "appears_codon_aware" in summary


def test_per_sequence_and_site_stats():
    profiler = MSATK(DATA / "dna.fasta")
    seq_rows = _rows(profiler.per_sequence_stats())
    site_rows = _rows(profiler.per_site_stats())
    assert len(seq_rows) == 3
    assert len(site_rows) == 9
    assert any(row["gap_count"] == 1 for row in seq_rows)
    assert any(row["gap_fraction"] > 0 for row in site_rows)


def test_codon_outputs():
    profiler = CodonProfiler(DATA / "dna.fasta")
    usage = profiler.codon_usage()
    rscu = profiler.rscu()
    stops = profiler.stop_codon_report()
    assert any(row["codon"] == "ATG" and row["count"] == 3 for row in usage)
    assert rscu
    assert any(row["codon"] == "TAA" for row in stops)


def test_protein_outputs():
    profiler = ProteinProfiler(DATA / "protein.faa")
    composition = profiler.amino_acid_composition()
    classes = profiler.residue_class_summary()
    hydropathy = profiler.hydrophobicity_summary()
    assert any(row["amino_acid"] == "M" and row["count"] == 3 for row in composition)
    assert any(row["residue_class"] == "hydrophobic" for row in classes)
    assert len(hydropathy) == 3


def test_write_outputs(tmp_path):
    profiler = MSATK(DATA / "dna.fasta")
    result = profiler.write_outputs(tmp_path, plots=False, force=True)
    assert result["summary"]["tool"] == "MSATK"
    assert (tmp_path / "tables" / "alignment_summary.csv").exists()
    assert (tmp_path / "tables" / "per_sequence_stats.csv").exists()
    assert (tmp_path / "report.html").exists()
    assert (tmp_path / "qc_warnings.txt").exists()
    assert (tmp_path / "summary.json").exists()
    assert (tmp_path / "assets" / "msatk_logo.png").exists()


def _rows(value):
    if hasattr(value, "to_dict"):
        return value.to_dict(orient="records")
    return value


def test_alignment_profiler_alias():
    assert AlignmentProfiler is MSATK


def test_lowercase_msatk_alias():
    assert msatk is MSATK
