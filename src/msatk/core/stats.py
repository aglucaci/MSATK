"""Backward-compatible import location for MSATK statistics (see :mod:`msatk.stats`)."""

from msatk.stats import (  # noqa: F401
    alignment_summary,
    composition_summary,
    distance_matrix,
    gc_content,
    identity_matrix,
    pairwise_identity,
    per_sequence_stats,
    per_site_stats,
    shannon_entropy,
)
