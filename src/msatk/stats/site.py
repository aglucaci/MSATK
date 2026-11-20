"""Per-site (column) statistics."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from msatk.core.alphabet import is_nucleotide
from msatk.core.encoded import (
    STATE_AMBIGUOUS,
    STATE_MISSING,
    EncodedAlignment,
)
from msatk.stats.information import entropy_from_counts


def _f(value: float) -> float | None:
    return (
        None if value is None or (isinstance(value, float) and math.isnan(value)) else float(value)
    )


def site_arrays(enc: EncodedAlignment) -> dict[str, np.ndarray]:
    """Vectorised per-site statistics (NaN marks undefined values)."""

    n = enc.n
    counts = enc.site_counts
    observed = counts.sum(axis=1)
    sorted_counts = np.sort(counts, axis=1)[:, ::-1]
    dominant = sorted_counts[:, 0] if enc.k else np.zeros(enc.width, dtype=int)
    states = (counts > 0).sum(axis=1)
    states_ge2 = (counts >= 2).sum(axis=1)
    entropy = entropy_from_counts(counts)
    with np.errstate(divide="ignore", invalid="ignore"):
        dom_freq = np.where(observed > 0, dominant / np.maximum(observed, 1), np.nan)
        norm_entropy = entropy / math.log2(enc.k) if enc.k > 1 else entropy * np.nan
    unobserved = observed == 0
    return {
        "observed": observed,
        "gaps": enc.gap_mask.sum(axis=0),
        "missing": (enc.state == STATE_MISSING).sum(axis=0),
        "ambiguous": (enc.state == STATE_AMBIGUOUS).sum(axis=0),
        "states": states,
        "dominant_count": dominant,
        "dominant_frequency": dom_freq,
        "consensus_support": dominant / n if n else dominant * 0.0,
        "entropy": entropy,
        "normalized_entropy": norm_entropy,
        "effective_states": np.where(unobserved, np.nan, np.power(2.0, entropy)),
        "information_content": np.where(unobserved, np.nan, math.log2(max(enc.k, 1)) - entropy),
        "unobserved": unobserved,
        "variable": states > 1,
        "constant": states == 1,
        "parsimony_informative": states_ge2 >= 2,
        "singleton": (states > 1) & (states_ge2 < 2),
    }


def per_site_stats_encoded(enc: EncodedAlignment) -> list[dict[str, Any]]:
    arrays = site_arrays(enc)
    n = enc.n
    alphabet = enc.policy.canonical
    consensus = enc.consensus_index
    codon = enc.policy.molecule_type == "codon"
    nucleotide = is_nucleotide(enc.policy.molecule_type)
    counts = enc.site_counts
    rows: list[dict[str, Any]] = []
    for i in range(enc.width):
        unobserved = bool(arrays["unobserved"][i])
        row: dict[str, Any] = {
            "site_index": i + 1,
            "site_class": "unobserved"
            if unobserved
            else ("variable" if arrays["variable"][i] else "constant"),
            "observed_count": int(arrays["observed"][i]),
            "occupancy": float(arrays["observed"][i] / n) if n else 0.0,
            "gap_count": int(arrays["gaps"][i]),
            "gap_fraction": float(arrays["gaps"][i] / n) if n else 0.0,
            "missing_fraction": float(arrays["missing"][i] / n) if n else 0.0,
            "ambiguous_fraction": float(arrays["ambiguous"][i] / n) if n else 0.0,
            "dominant_character": alphabet[consensus[i]] if consensus[i] >= 0 else None,
            "dominant_character_frequency": _f(arrays["dominant_frequency"][i]),
            "conservation_score": _f(arrays["dominant_frequency"][i]),
            "consensus_support": float(arrays["consensus_support"][i]),
            "observed_state_count": int(arrays["states"][i]),
            "entropy": _f(arrays["entropy"][i]),
            "normalized_entropy": _f(arrays["normalized_entropy"][i]),
            "effective_states": _f(arrays["effective_states"][i]),
            "information_content": _f(arrays["information_content"][i]),
            "constant": bool(arrays["constant"][i]),
            "variable": bool(arrays["variable"][i]),
            "parsimony_informative": bool(arrays["parsimony_informative"][i]),
            "singleton": bool(arrays["singleton"][i]),
        }
        if nucleotide:
            observed = counts[i].sum()
            gc = counts[i][alphabet.index("G")] + counts[i][alphabet.index("C")]
            row["gc_content"] = float(gc / observed) if observed else None
        if codon:
            row["codon_position"] = (i % 3) + 1
        rows.append(row)
    return rows
