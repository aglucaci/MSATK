"""Substitution-matrix similarity for protein alignments (BLOSUM62).

The matrix is bundled as ``msatk/data/blosum62.json`` (NCBI BLOSUM62).
"""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from typing import Any

import numpy as np

from msatk.core.encoded import EncodedAlignment


@lru_cache(maxsize=1)
def blosum62() -> tuple[str, np.ndarray]:
    data = json.loads(resources.files("msatk.data").joinpath("blosum62.json").read_text())
    return str(data["alphabet"]), np.array(data["scores"], dtype=np.int16)


def blosum62_score(a: str, b: str) -> int | None:
    """BLOSUM62 score for two residues, or ``None`` when either is not in the matrix."""

    alphabet, matrix = blosum62()
    i, j = alphabet.find(a.upper()), alphabet.find(b.upper())
    return None if i < 0 or j < 0 else int(matrix[i, j])


def _canonical_matrix(canonical: str) -> np.ndarray:
    alphabet, matrix = blosum62()
    idx = [alphabet.index(ch) for ch in canonical]
    return matrix[np.ix_(idx, idx)].astype(np.float64)


def column_substitution_scores(enc: EncodedAlignment) -> list[dict[str, Any]]:
    """Mean BLOSUM62 score over all pairs of canonical residues in each column.

    With residue counts ``c`` the sum over unordered pairs of distinct sequences is
    ``(cᵀ S c − Σ_a c_a S_aa) / 2`` and the number of such pairs is ``n(n−1)/2``.
    Columns with fewer than two residues are undefined (``None``).
    """

    s = _canonical_matrix(enc.policy.canonical)
    counts = enc.site_counts.astype(np.float64)
    n = counts.sum(axis=1)
    total = np.einsum("la,ab,lb->l", counts, s, counts) - counts @ np.diag(s)
    pairs = n * (n - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(pairs > 0, total / np.where(pairs > 0, pairs, 1), np.nan)
    return [
        {
            "site_index": i + 1,
            "residue_pairs": int(p // 2),
            "mean_blosum62": None if np.isnan(m) else float(m),
        }
        for i, (p, m) in enumerate(zip(pairs, mean))
    ]


def substitution_summary(enc: EncodedAlignment) -> dict[str, Any]:
    rows = column_substitution_scores(enc)
    values = [r["mean_blosum62"] for r in rows if r["mean_blosum62"] is not None]
    return {
        "mean_column_blosum62": float(np.mean(values)) if values else None,
        "median_column_blosum62": float(np.median(values)) if values else None,
        "columns_scored_blosum62": len(values),
    }
