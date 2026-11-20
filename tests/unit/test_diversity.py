"""Population-genetic statistics validated against hand calculation and scikit-allel."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pytest

from msatk.models import Alignment, SequenceRecord
from msatk.stats import compute_statistics, tajimas_d


def _aln(seqs):
    return Alignment(tuple(SequenceRecord(f"s{i}", s) for i, s in enumerate(seqs)))


def _random_population(seed=3, n=12, length=200, mutations=25):
    rng = np.random.default_rng(seed)
    anc = rng.integers(0, 4, length)
    m = np.tile(anc, (n, 1))
    for _ in range(mutations):
        site = rng.integers(0, length)
        carriers = rng.random(n) < rng.random()
        m[carriers, site] = (anc[site] + rng.integers(1, 4)) % 4
    return ["".join("ACGT"[x] for x in row) for row in m]


def test_pi_and_theta_match_brute_force():
    seqs = _random_population()
    s = compute_statistics(_aln(seqs), "dna").summary
    n, length = len(seqs), len(seqs[0])
    diffs = [sum(a != b for a, b in zip(x, y)) for x, y in combinations(seqs, 2)]
    pi = sum(diffs) / len(diffs) / length
    seg = sum(len({x[i] for x in seqs}) > 1 for i in range(length))
    a1 = sum(1 / i for i in range(1, n))
    assert s["nucleotide_diversity_pi"] == pytest.approx(pi)
    assert s["segregating_sites"] == seg
    assert s["theta_watterson"] == pytest.approx(seg / a1 / length)


def test_tajimas_d_reference_values():
    # Values cross-checked against scikit-allel 1.3 (allel.tajima_d) for these alignments.
    s = compute_statistics(_aln(_random_population(seed=3)), "dna").summary
    assert s["tajimas_d"] == pytest.approx(0.75264128, abs=1e-7)
    assert tajimas_d(3, 5, 2.0) is None  # n < 4
    assert tajimas_d(10, 0, 0.0) is None  # no segregating sites


def test_against_scikit_allel_when_available():
    allel = pytest.importorskip("allel")
    for seed in (1, 2, 5):
        seqs = _random_population(seed=seed)
        m = np.array([["ACGT".index(c) for c in s] for s in seqs])
        ac = allel.AlleleCountsArray(np.stack([(m == c).sum(0) for c in range(4)], 1))
        pos = np.arange(1, m.shape[1] + 1)
        s = compute_statistics(_aln(seqs), "dna").summary
        assert s["nucleotide_diversity_pi"] == pytest.approx(
            allel.sequence_diversity(pos, ac, start=1, stop=m.shape[1])
        )
        assert s["tajimas_d"] == pytest.approx(allel.tajima_d(ac))


def test_haplotypes_and_tstv():
    # Two haplotypes (3x AAAA, 1x GACA): one transition (A<->G) and one transversion (A<->C)
    seqs = ["AAAA", "AAAA", "AAAA", "GACA"]
    s = compute_statistics(_aln(seqs), "dna").summary
    assert s["haplotype_count"] == 2
    assert s["haplotype_diversity"] == pytest.approx(4 / 3 * (1 - (0.75**2 + 0.25**2)))
    assert s["transitions_pairwise"] == 3
    assert s["transversions_pairwise"] == 3
    assert s["ts_tv_ratio"] == pytest.approx(1.0)


def test_popgen_uses_complete_sites_only():
    s = compute_statistics(_aln(["AC-T", "ACGT", "ATGT"]), "dna").summary
    assert s["popgen_sites_used"] == 3
    assert s["segregating_sites"] == 1


def test_protein_has_no_popgen_fields():
    s = compute_statistics(_aln(["MKV", "MKI"]), "protein").summary
    assert "nucleotide_diversity_pi" not in s
    assert s["gc_content"] is None
