"""Population-genetic and evolutionary diversity statistics for nucleotide alignments.

These statistics assume the sequences are homologous samples from a population
(or at least a set where the infinite-sites interpretation is sensible). MSATK
computes them on *complete sites* only — columns in which every sequence carries
a canonical base — and reports how many sites were used. They are reported for
nucleotide alignments only and should be interpreted with care for arbitrary
multi-species alignments.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from msatk.core.alphabet import is_nucleotide
from msatk.core.encoded import EncodedAlignment

# Canonical order is ACGT (or ACGU): purines A(0), G(2); pyrimidines C(1), T/U(3).
_TRANSITION_PAIRS = ((0, 2), (1, 3))
_TRANSVERSION_PAIRS = ((0, 1), (0, 3), (2, 1), (2, 3))


def tajimas_d(n: int, segregating: int, pi_total: float) -> float | None:
    """Tajima's (1989) D from sample size, segregating sites and total mean pairwise differences."""

    if n < 4 or segregating == 0:
        return None
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))
    b1 = (n + 1) / (3.0 * (n - 1))
    b2 = 2.0 * (n * n + n + 3) / (9.0 * n * (n - 1))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    denom = math.sqrt(e1 * segregating + e2 * segregating * (segregating - 1))
    return (pi_total - segregating / a1) / denom if denom else None


def diversity_summary(enc: EncodedAlignment) -> dict[str, Any]:
    """Nucleotide diversity, Watterson's theta, Tajima's D, haplotype diversity and Ts/Tv."""

    if not is_nucleotide(enc.policy.molecule_type):
        return {}
    n = enc.n
    complete = enc.residue_mask.all(axis=0)
    counts = enc.site_counts[complete]
    used = int(complete.sum())
    result: dict[str, Any] = {
        "popgen_sites_used": used,
        "segregating_sites": None,
        "nucleotide_diversity_pi": None,
        "theta_watterson": None,
        "tajimas_d": None,
        "haplotype_count": None,
        "haplotype_diversity": None,
        "transitions_pairwise": None,
        "transversions_pairwise": None,
        "ts_tv_ratio": None,
    }
    if n < 2 or used == 0:
        return result
    segregating = int(((counts > 0).sum(axis=1) > 1).sum())
    pairs = n * (n - 1) / 2.0
    # Number of differing pairs at each site = C(n,2) - sum_c C(c,2)
    same_pairs = (counts * (counts - 1) / 2.0).sum(axis=1)
    diff_pairs = pairs - same_pairs
    pi_total = float(diff_pairs.sum() / pairs)
    a1 = sum(1.0 / i for i in range(1, n))
    ts = sum(float((counts[:, a] * counts[:, b]).sum()) for a, b in _TRANSITION_PAIRS)
    tv = sum(float((counts[:, a] * counts[:, b]).sum()) for a, b in _TRANSVERSION_PAIRS)
    haplotypes: dict[bytes, int] = {}
    sub = enc.chars[:, complete]
    for row in sub:
        key = row.tobytes()
        haplotypes[key] = haplotypes.get(key, 0) + 1
    freqs = np.array(list(haplotypes.values()), dtype=float) / n
    result.update(
        {
            "segregating_sites": segregating,
            "nucleotide_diversity_pi": pi_total / used,
            "theta_watterson": segregating / a1 / used,
            "tajimas_d": tajimas_d(n, segregating, pi_total),
            "haplotype_count": len(haplotypes),
            "haplotype_diversity": float(n / (n - 1) * (1.0 - (freqs**2).sum())),
            "transitions_pairwise": int(ts),
            "transversions_pairwise": int(tv),
            "ts_tv_ratio": ts / tv if tv else None,
        }
    )
    return result
