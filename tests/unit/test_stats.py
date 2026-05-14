import pytest

from msatk.core.stats import gc_content, per_sequence_stats, per_site_stats, shannon_entropy
from msatk.io import read_alignment


def test_entropy_conserved_site():
    assert shannon_entropy(["A", "A", "A"]) == 0.0


def test_entropy_variable_site():
    assert shannon_entropy(["A", "T"]) == pytest.approx(1.0)


def test_gc_content():
    assert gc_content("ATGC") == 0.5


def test_gap_statistics(data_dir):
    alignment = read_alignment(data_dir / "dna" / "dna_with_gaps.fasta")
    seq_rows = per_sequence_stats(alignment)
    site_rows = per_site_stats(alignment)
    assert any(row["gap_fraction"] > 0 for row in seq_rows)
    assert any(row["gap_fraction"] > 0 for row in site_rows)
