"""Column diversity beyond Shannon entropy, and site-pattern statistics.

* Sequence-weighted entropy down-weights clusters of near-identical sequences
  (weights from :mod:`msatk.stats.redundancy`).
* Simpson concentration ``sum p^2`` and Gini–Simpson ``1 - sum p^2`` are
  alternative column-diversity measures that are less sensitive to rare states.
* Site patterns are whole columns. Likelihood software compresses identical
  columns; the number of distinct patterns and their multiplicities describe how
  much distinct character information the alignment holds.
* Variable sites are classed by the number of observed states and by the count of
  the minor (non-majority) states: singleton (1), doubleton (2) and so on.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from msatk.core.encoded import EncodedAlignment


def _entropy_bits(p: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0.0)
    return -(p * logs).sum(axis=1)


def weighted_site_counts(enc: EncodedAlignment, weights: np.ndarray) -> np.ndarray:
    """``(L, k)`` weighted residue counts."""

    out = np.zeros((enc.width, enc.k), dtype=np.float64)
    for c in range(enc.k):
        out[:, c] = weights @ (enc.residue_index == c)
    return out


def site_pattern_analysis(
    enc: EncodedAlignment, weights: np.ndarray | None = None
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    counts = enc.site_counts.astype(np.float64)
    observed = counts.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(observed[:, None] > 0, counts / np.maximum(observed[:, None], 1), 0.0)
    simpson = (p**2).sum(axis=1)
    entropy = _entropy_bits(p)
    unobserved = observed == 0
    w_entropy = None
    if weights is not None:
        wc = weighted_site_counts(enc, weights)
        wt = wc.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            wp = np.where(wt[:, None] > 0, wc / np.maximum(wt[:, None], 1e-300), 0.0)
        w_entropy = _entropy_bits(wp)

    # Site patterns: identical full columns (normalised characters, gaps included).
    cols = np.ascontiguousarray(enc.chars.T)
    keys = [cols[j].tobytes() for j in range(enc.width)]
    multiplicity: dict[bytes, int] = {}
    first_site: dict[bytes, int] = {}
    for j, key in enumerate(keys):
        multiplicity[key] = multiplicity.get(key, 0) + 1
        first_site.setdefault(key, j)
    mult = np.array(list(multiplicity.values()), dtype=np.float64)
    pattern_p = mult / mult.sum() if mult.size else mult
    pattern_entropy = float(-(pattern_p * np.log2(pattern_p)).sum()) if mult.size else 0.0

    states = (counts > 0).sum(axis=1)
    minor = observed - counts.max(axis=1, initial=0)
    variable = states >= 2
    state_hist = {
        "sites_1_state": int((states == 1).sum()),
        "sites_2_states": int((states == 2).sum()),
        "sites_3_states": int((states == 3).sum()),
        "sites_4plus_states": int((states >= 4).sum()),
    }
    minor_hist = {
        "singleton_minor_sites": int((variable & (minor == 1)).sum()),
        "doubleton_minor_sites": int((variable & (minor == 2)).sum()),
        "minor_count_3plus_sites": int((variable & (minor >= 3)).sum()),
    }
    obs_sites = ~unobserved
    summary: dict[str, Any] = {
        "distinct_site_patterns": len(multiplicity),
        "site_pattern_ratio": len(multiplicity) / enc.width if enc.width else None,
        "repeated_pattern_sites": int(enc.width - len(multiplicity)),
        "max_pattern_multiplicity": int(mult.max()) if mult.size else 0,
        "site_pattern_entropy": pattern_entropy,
        "mean_gini_simpson": float((1 - simpson)[obs_sites].mean()) if obs_sites.any() else None,
        "mean_effective_alphabet_size": float(np.power(2.0, entropy)[obs_sites].mean())
        if obs_sites.any()
        else None,
        "mean_weighted_entropy": float(w_entropy[obs_sites].mean())
        if w_entropy is not None and obs_sites.any()
        else None,
        "weighted": weights is not None,
        **state_hist,
        **minor_hist,
    }
    rows = []
    for j in range(enc.width):
        unobs = bool(unobserved[j])
        rows.append(
            {
                "site_index": j + 1,
                "simpson": None if unobs else float(simpson[j]),
                "gini_simpson": None if unobs else float(1 - simpson[j]),
                "effective_alphabet_size": None if unobs else float(2 ** entropy[j]),
                "weighted_entropy": None if unobs or w_entropy is None else float(w_entropy[j]),
                "minor_state_count": int(minor[j]),
                "pattern_multiplicity": multiplicity[keys[j]],
                "pattern_first_site": first_site[keys[j]] + 1,
            }
        )
    return summary, {"site_diversity": rows}


def jensen_shannon(p: np.ndarray, q: np.ndarray) -> float:
    """Jensen–Shannon divergence in bits (0..1)."""

    m = 0.5 * (p + q)

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        mask = a > 0
        return float((a[mask] * np.log2(a[mask] / b[mask])).sum())

    return 0.5 * kl(p, m) + 0.5 * kl(q, m) if math.isfinite(m.sum()) else float("nan")
