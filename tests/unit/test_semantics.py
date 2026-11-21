"""Hand-checkable tests of MSATK's statistical semantics (gaps, missing data, undefined values)."""

from __future__ import annotations

import pytest

from msatk.core.alphabet import policy_for
from msatk.models import Alignment, SequenceRecord
from msatk.stats import compute_statistics, pairwise_identity, shannon_entropy


def aln(*seqs: str) -> Alignment:
    return Alignment(tuple(SequenceRecord(f"s{i + 1}", s) for i, s in enumerate(seqs)))


def test_all_gap_column_is_unobserved_not_conserved(data_dir):
    from msatk.io import read_alignment

    bundle = compute_statistics(
        read_alignment(data_dir / "edge_cases" / "all_gap_column.fasta"), "dna"
    )
    site4 = bundle.sites[3]
    assert site4["site_class"] == "unobserved"
    assert site4["constant"] is False and site4["variable"] is False
    assert site4["entropy"] is None and site4["conservation_score"] is None
    assert bundle.summary["unobserved_sites"] == 1
    assert bundle.summary["constant_sites"] == 6
    assert bundle.summary["conserved_sites"] == 6


def test_identity_undefined_without_comparable_sites():
    bundle = compute_statistics(aln("AC--", "--GT", "ACGT"), "dna")
    pw = bundle.pairwise
    assert pw is not None and pw.is_exact
    assert pw.identity[0, 1] != pw.identity[0, 1]  # NaN
    assert pw.pairs_without_comparable_sites == 1
    assert bundle.summary["mean_pairwise_identity"] == pytest.approx(1.0)
    assert pairwise_identity("AC--", "--GT") is None
    long = {(r["sequence_a"], r["sequence_b"]): r for r in pw.long_table()}
    assert long[("s1", "s2")]["identity"] is None
    assert long[("s1", "s2")]["comparable_sites"] == 0


def test_padding_is_missing_not_gap(data_dir):
    from msatk.io import read_alignment

    a = read_alignment(data_dir / "dna" / "unequal_lengths.fasta")
    bundle = compute_statistics(a, "dna")
    assert bundle.summary["padded_cells"] == 4
    assert bundle.summary["total_gaps"] == 0
    seq2 = bundle.sequences[1]
    assert seq2["padded_cells"] == 4 and seq2["gap_count"] == 0
    assert seq2["missing_count"] == 0  # padding reported separately
    assert bundle.sites[-1]["occupancy"] == pytest.approx(0.5)


def test_strict_mode_rejects_unequal_lengths(data_dir):
    from msatk.exceptions import AlignmentValidationError
    from msatk.io import read_alignment

    with pytest.raises(AlignmentValidationError):
        read_alignment(data_dir / "dna" / "unequal_lengths.fasta", validation_mode="strict")


def test_asparagine_is_a_residue_in_proteins_but_N_is_missing_in_dna():
    assert policy_for("protein").state("N") == "residue"
    assert policy_for("dna").state("N") == "missing"
    assert policy_for("protein").state("X") == "missing"
    assert policy_for("dna").state("R") == "ambiguous"
    prot = compute_statistics(aln("NNNN", "NNNN"), "protein")
    assert prot.summary["constant_sites"] == 4
    assert prot.summary["residue_cells"] == 8
    dna = compute_statistics(aln("NNNN", "NNNN"), "dna")
    assert dna.summary["unobserved_sites"] == 4


def test_entropy_definitions():
    assert shannon_entropy(["A", "A", "A"]) == 0.0
    assert shannon_entropy(["A", "T"]) == pytest.approx(1.0)
    assert shannon_entropy(["-", "-"], "dna") is None
    # Ambiguity codes and gaps are excluded from the distribution.
    assert shannon_entropy(["A", "R", "-", "N"], "dna") == 0.0
    bundle = compute_statistics(aln("A", "C", "G", "T"), "dna")
    site = bundle.sites[0]
    assert site["entropy"] == pytest.approx(2.0)
    assert site["normalized_entropy"] == pytest.approx(1.0)
    assert site["effective_states"] == pytest.approx(4.0)
    assert site["information_content"] == pytest.approx(0.0)


def test_parsimony_and_singleton_definitions():
    # col1 constant; col2 singleton (one T); col3 PI (AA/TT); col4 PI with extra singleton (AA/TT/G)
    bundle = compute_statistics(aln("AAAA", "AAAA", "ATTT", "AATT", "AATG"), "dna")
    classes = [(s["parsimony_informative"], s["singleton"]) for s in bundle.sites]
    assert classes == [(False, False), (False, True), (True, False), (True, False)]
    assert bundle.summary["parsimony_informative_sites"] == 2
    assert bundle.summary["singleton_sites"] == 1


def test_gap_runs_and_spans():
    bundle = compute_statistics(aln("--AC-GT---", "ACACAGTACA"), "dna")
    row = bundle.sequences[0]
    assert row["leading_gaps"] == 2
    assert row["trailing_gaps"] == 3
    assert row["internal_gap_count"] == 1
    assert row["gap_run_count"] == 3
    assert row["max_gap_run"] == 3
    assert (row["span_start"], row["span_end"], row["span_length"]) == (3, 7, 5)
    assert row["ungapped_length"] == 4


def test_gc_pooled_and_skew():
    bundle = compute_statistics(aln("GGGC", "AAAT"), "dna")
    assert bundle.summary["gc_content"] == pytest.approx(0.5)
    assert bundle.summary["gc_skew"] == pytest.approx(0.5)
    assert bundle.summary["at_skew"] == pytest.approx(0.5)
    assert bundle.sequences[0]["gc_content"] == pytest.approx(1.0)


def test_codon_position_only_for_codon_alignments():
    dna = compute_statistics(aln("ATGAAA", "ATGAAG"), "dna")
    codon = compute_statistics(aln("ATGAAA", "ATGAAG"), "codon")
    assert "codon_position" not in dna.sites[0]
    assert [s["codon_position"] for s in codon.sites] == [1, 2, 3, 1, 2, 3]
    prot = compute_statistics(aln("MKV", "MKI"), "protein")
    assert "codon_position" not in prot.sites[0]
    assert "gc_content" not in prot.sites[0]


def test_duplicates_and_consensus_identity():
    bundle = compute_statistics(aln("ACGT", "ACGT", "ACGA"), "dna")
    assert bundle.sequences[1]["duplicate_of"] == "s1"
    assert bundle.summary["duplicate_sequence_count"] == 1
    assert bundle.summary["unique_sequence_count"] == 2
    assert bundle.sequences[2]["consensus_identity"] == pytest.approx(0.75)


def test_sampled_pairwise_is_deterministic_and_close_to_exact():
    import numpy as np

    rng = np.random.default_rng(0)
    seqs = ["".join(rng.choice(list("ACGT"), 50)) for _ in range(40)]
    exact = compute_statistics(aln(*seqs), "dna", pairwise="exact").summary
    s1 = compute_statistics(aln(*seqs), "dna", pairwise="sampled", sample_pairs=300, seed=3).summary
    s2 = compute_statistics(aln(*seqs), "dna", pairwise="sampled", sample_pairs=300, seed=3).summary
    assert s1["pairwise_method"] == "sampled"
    assert s1["mean_pairwise_identity"] == s2["mean_pairwise_identity"]
    assert abs(s1["mean_pairwise_identity"] - exact["mean_pairwise_identity"]) < 0.02
    none = compute_statistics(aln(*seqs), "dna", pairwise="none").summary
    assert none["mean_pairwise_identity"] is None and none["pairwise_method"] == "none"


def test_auto_pairwise_switches_to_sampling_above_threshold():
    seqs = ["ACGT"] * 12
    summary = compute_statistics(aln(*seqs), "dna", pairwise="auto", max_exact_sequences=10).summary
    assert summary["pairwise_method"] in {"sampled", "all_pairs"}


def test_single_sequence_has_undefined_pairwise():
    bundle = compute_statistics(aln("ACGT"), "dna")
    assert bundle.pairwise is None
    assert bundle.summary["mean_pairwise_identity"] is None
    assert bundle.summary["nucleotide_diversity_pi"] is None
