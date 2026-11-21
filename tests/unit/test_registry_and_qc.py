"""Every emitted field must be defined in the metric dictionary; QC rules fire as documented."""

from __future__ import annotations

import pytest

from msatk import MSATK
from msatk.core.qc import QCThresholds, qc_findings
from msatk.models import Alignment, SequenceRecord
from msatk.stats import compute_statistics
from msatk.stats.registry import is_defined, metric_dictionary, metrics_markdown

FIXTURES = [
    ("dna/tiny_dna.fasta", "dna"),
    ("dna/dna_with_gaps.fasta", "dna"),
    ("codon/valid_codon_alignment.fasta", "codon"),
    ("protein/protein_with_gaps.faa", "protein"),
    ("dna/unequal_lengths.fasta", "dna"),
]


@pytest.mark.parametrize("path,kind", FIXTURES)
def test_all_fields_defined(data_dir, path, kind):
    profiler = MSATK(data_dir / path, sequence_type=kind)
    bundle = profiler.stats()
    summary = profiler._provenance_summary()
    missing = [k for k in summary if not is_defined(k, "alignment")]
    missing += [k for r in bundle.sites for k in r if not is_defined(k, "site")]
    missing += [k for r in bundle.sequences for k in r if not is_defined(k, "sequence")]
    if bundle.pairwise is not None:
        missing += [k for r in bundle.pairwise.long_table() for k in r if not is_defined(k, "pair")]
    assert not sorted(set(missing))


def test_dictionary_exports():
    entries = metric_dictionary()
    assert len(entries) > 100
    assert all({"name", "level", "description"} <= set(e) for e in entries)
    assert "| `nucleotide_diversity_pi` |" in metrics_markdown()


def _aln(seqs, ids=None):
    ids = ids or [f"s{i}" for i in range(len(seqs))]
    return Alignment(tuple(SequenceRecord(i, s) for i, s in zip(ids, seqs)))


def _codes(a, kind, **kw):
    b = compute_statistics(a, kind)
    return {
        f.code: f for f in qc_findings(a, kind, b.sequences, b.sites, b.summary, QCThresholds(**kw))
    }


def test_qc_rules():
    clean = _codes(_aln(["ACGTACGTAC"] * 3 + ["ACGTACGTAA"]), "dna")
    assert "NO_MAJOR_ISSUES" in clean
    codes = _codes(
        _aln(
            ["ACGT----AC", "ACGTACGTAC", "TTTTTTTTTT", "ACGTACGTAC", "-CGTACGTAC"],
            ids=["a", "b", "c", "b", "e"],
        ),
        "dna",
    )
    assert {
        "DUPLICATE_IDS",
        "DUPLICATE_SEQUENCES",
        "DIVERGENT_SEQUENCES",
        "HIGH_GAP_SEQUENCES",
    } <= set(codes)
    assert codes["DIVERGENT_SEQUENCES"].items == ["c"]
    gaps = _codes(_aln(["A-", "A-", "AC"]), "dna", max_gap_site=0.5)
    assert gaps["HIGH_GAP_SITES"].items == ["2"]
    nonoverlap = _codes(_aln(["AC--", "--GT"]), "dna")
    assert "NON_OVERLAPPING_PAIRS" in nonoverlap
    codon = _codes(_aln(["ATGTAAGGG", "ATG-AAGGG"]), "codon")
    assert {"INTERNAL_STOP_CODON", "PARTIAL_CODON_GAPS"} <= set(codon)
    prot = _codes(_aln(["MK*V", "MKAV"]), "protein")
    assert "INTERNAL_STOPS_PROTEIN" in prot


def test_findings_sorted_by_severity():
    a = _aln(["ACGT", "AC--", "ACGT"], ids=["x", "y", "x"])
    b = compute_statistics(a, "dna")
    sev = [f.severity for f in qc_findings(a, "dna", b.sequences, b.sites, b.summary)]
    order = {"error": 0, "warning": 1, "info": 2, "ok": 3}
    assert sev == sorted(sev, key=order.get)


def test_docs_metrics_page_is_current():
    from pathlib import Path

    page = Path(__file__).resolve().parents[2] / "docs" / "metrics.md"
    assert page.read_text().strip() == metrics_markdown().strip(), (
        "Regenerate with: msatk metrics > docs/metrics.md"
    )
