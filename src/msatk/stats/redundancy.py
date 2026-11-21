"""Sequence redundancy: effective number of sequences, weights and clusters.

An alignment of 10,000 near-identical sequences carries far less independent
information than its row count suggests. The effective number of sequences is

    N_eff = sum_i 1 / n_i,

where ``n_i`` is the number of sequences (including ``i``) whose identity to ``i``
is at least ``threshold``; ``w_i = 1 / n_i`` is the sequence weight used by
coevolution and profile methods (e.g. Morcos et al. 2011; Hopf et al. 2017).

Identity is computed over sites where both sequences carry a canonical residue.
A pair counts as neighbours only when its comparable sites cover at least
``min_overlap`` of the *longer* sequence's residues, so a short fragment is not
counted as identical to every sequence it happens to match. Pairs without enough
overlap are treated as non-neighbours.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from msatk.core.encoded import EncodedAlignment
from msatk.stats.pairwise import PairwiseResult

DEFAULT_THRESHOLD = 0.8
DEFAULT_MIN_OVERLAP = 0.5


def neighbour_matrix(
    enc: EncodedAlignment,
    pw: PairwiseResult,
    threshold: float = DEFAULT_THRESHOLD,
    min_overlap: float = DEFAULT_MIN_OVERLAP,
) -> np.ndarray:
    """Boolean ``(N, N)`` matrix of neighbours (diagonal True)."""

    if pw.identity is None or pw.comparable is None:
        raise ValueError("Sequence weights need exact pairwise identities.")
    residues = enc.residue_mask.sum(axis=1).astype(np.float64)
    longer = np.maximum(residues[:, None], residues[None, :])
    overlap_ok = pw.comparable >= np.maximum(1.0, min_overlap * longer)
    ident = np.nan_to_num(pw.identity, nan=-1.0)
    nb = (ident >= threshold) & overlap_ok
    np.fill_diagonal(nb, True)
    return nb


def sequence_weights(
    enc: EncodedAlignment,
    pw: PairwiseResult,
    threshold: float = DEFAULT_THRESHOLD,
    min_overlap: float = DEFAULT_MIN_OVERLAP,
) -> np.ndarray:
    """``w_i = 1 / n_i``; sums to N_eff."""

    return 1.0 / neighbour_matrix(enc, pw, threshold, min_overlap).sum(axis=1)


def _components(nb: np.ndarray) -> np.ndarray:
    """Connected components (single linkage at the threshold); labels in first-seen order."""

    n = nb.shape[0]
    labels = np.full(n, -1, dtype=np.int64)
    current = 0
    for start in range(n):
        if labels[start] >= 0:
            continue
        stack = [start]
        labels[start] = current
        while stack:
            i = stack.pop()
            for j in np.flatnonzero(nb[i] & (labels < 0)):
                labels[j] = current
                stack.append(int(j))
        current += 1
    return labels


def redundancy_analysis(
    enc: EncodedAlignment,
    pw: PairwiseResult,
    threshold: float = DEFAULT_THRESHOLD,
    min_overlap: float = DEFAULT_MIN_OVERLAP,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]], np.ndarray]:
    """Summary, per-sequence weights/clusters, and the weight vector."""

    nb = neighbour_matrix(enc, pw, threshold, min_overlap)
    n_i = nb.sum(axis=1)
    weights = 1.0 / n_i
    neff = float(weights.sum())
    n = enc.n
    labels = _components(nb)
    sizes = np.bincount(labels)
    identities = pw.pair_identities[~np.isnan(pw.pair_identities)]
    qs = (0.05, 0.25, 0.5, 0.75, 0.95)
    summary: dict[str, Any] = {
        "identity_threshold": threshold,
        "min_overlap": min_overlap,
        "effective_sequences": neff,
        "effective_sequence_fraction": neff / n if n else None,
        "effective_sequences_per_site": neff / enc.width if enc.width else None,
        "mean_neighbours": float(n_i.mean()),
        "max_neighbours": int(n_i.max()),
        "clusters": int(sizes.size),
        "largest_cluster": int(sizes.max()),
        "singleton_clusters": int((sizes == 1).sum()),
        "pairwise_identity_quantiles": {
            f"q{int(q * 100):02d}": float(np.quantile(identities, q)) if identities.size else None
            for q in qs
        },
    }
    rows = [
        {
            "sequence_id": seq_id,
            "neighbours": int(n_i[i]),
            "weight": float(weights[i]),
            "cluster": int(labels[i]),
            "cluster_size": int(sizes[labels[i]]),
        }
        for i, seq_id in enumerate(enc.alignment.ids)
    ]
    size_hist: dict[int, int] = {}
    for s in sizes.tolist():
        size_hist[s] = size_hist.get(s, 0) + 1
    tables: dict[str, list[dict[str, Any]]] = {
        "sequence_weights": rows,
        "cluster_size_distribution": [
            {"cluster_size": k, "clusters": v} for k, v in sorted(size_hist.items())
        ],
    }
    return summary, tables, weights
