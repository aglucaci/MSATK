"""Alignment statistics at alignment, sequence, site and pair resolution.

All functions share one vectorised encoding (:mod:`msatk.core.encoded`) and one
character policy (:mod:`msatk.core.alphabet`), so gaps, missing data and
ambiguity codes are treated identically by every metric. Every reported field is
defined in :mod:`msatk.stats.registry`.
"""

from __future__ import annotations

from typing import Any

from msatk.core.alphabet import DNA_POLICY
from msatk.core.detect import detect_molecule_type
from msatk.core.encoded import encode_alignment
from msatk.models import Alignment
from msatk.stats.alignment import (
    StatisticsBundle,
    composition_rows,
    compute_statistics,
    gap_block_rows,
)
from msatk.stats.diversity import diversity_summary, tajimas_d
from msatk.stats.information import entropy_from_counts, shannon_entropy
from msatk.stats.pairwise import PairwiseResult, compute_pairwise, pairwise_identity
from msatk.stats.sequence import per_sequence_stats_encoded
from msatk.stats.site import per_site_stats_encoded


def _type(alignment: Alignment, molecule_type: str | None) -> str:
    return molecule_type or detect_molecule_type(alignment)


def per_site_stats(alignment: Alignment, molecule_type: str | None = None) -> list[dict[str, Any]]:
    return per_site_stats_encoded(encode_alignment(alignment, _type(alignment, molecule_type)))


def per_sequence_stats(
    alignment: Alignment, molecule_type: str | None = None, pairwise: str = "auto"
) -> list[dict[str, Any]]:
    enc = encode_alignment(alignment, _type(alignment, molecule_type))
    return per_sequence_stats_encoded(enc, compute_pairwise(enc, method=pairwise))


def alignment_summary(
    alignment: Alignment, molecule_type: str | None = None, pairwise: str = "auto"
) -> dict[str, Any]:
    return compute_statistics(alignment, _type(alignment, molecule_type), pairwise=pairwise).summary


def composition_summary(
    alignment: Alignment, molecule_type: str | None = None
) -> list[dict[str, Any]]:
    return composition_rows(encode_alignment(alignment, _type(alignment, molecule_type)))


def identity_matrix(
    alignment: Alignment, molecule_type: str | None = None
) -> list[list[float | None]]:
    enc = encode_alignment(alignment, _type(alignment, molecule_type))
    result = compute_pairwise(enc, method="exact")
    if result is None or result.identity is None:
        return [[1.0]]
    import numpy as np

    np.fill_diagonal(result.identity, 1.0)
    return result.matrix_rows("identity")


def distance_matrix(
    alignment: Alignment, molecule_type: str | None = None
) -> list[list[float | None]]:
    return [
        [None if v is None else 1.0 - v for v in row]
        for row in identity_matrix(alignment, molecule_type)
    ]


def gc_content(sequence: str) -> float | None:
    """GC fraction among canonical bases (gaps, N and ambiguity codes excluded)."""

    bases = [DNA_POLICY.normalize(ch) for ch in sequence if DNA_POLICY.is_residue(ch)]
    if not bases:
        return None
    return (bases.count("G") + bases.count("C")) / len(bases)


__all__ = [
    "PairwiseResult",
    "StatisticsBundle",
    "alignment_summary",
    "composition_summary",
    "compute_pairwise",
    "compute_statistics",
    "distance_matrix",
    "diversity_summary",
    "entropy_from_counts",
    "gap_block_rows",
    "gc_content",
    "identity_matrix",
    "pairwise_identity",
    "per_sequence_stats",
    "per_site_stats",
    "shannon_entropy",
    "tajimas_d",
]
