from msatk.codon.analysis import codon_usage, iter_codons, rscu, stop_codon_report
from msatk.io import read_alignment


def test_split_codons():
    assert [codon for _, codon in iter_codons("ATGAAATTT")] == ["ATG", "AAA", "TTT"]


def test_codon_usage(data_dir):
    alignment = read_alignment(data_dir / "codon" / "valid_codon_alignment.fasta")
    usage = codon_usage(alignment)
    assert any(row["codon"] == "ATG" and row["count"] == 3 for row in usage)
    assert rscu(alignment)


def test_stop_codon_detection(data_dir):
    alignment = read_alignment(data_dir / "codon" / "codon_with_stop.fasta")
    stops = stop_codon_report(alignment)
    assert any(row["sequence_id"] == "seq1" for row in stops)
