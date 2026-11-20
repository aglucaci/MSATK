"""Information-theoretic site statistics."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable

import numpy as np

from msatk.core.alphabet import AlphabetPolicy, policy_for


def shannon_entropy(
    values: Iterable[str], molecule_type: str | AlphabetPolicy | None = None
) -> float | None:
    """Shannon entropy (bits) of the canonical residues in ``values``.

    Gaps, missing data and ambiguity codes are excluded from the distribution.
    Returns ``None`` when no canonical residue is observed: the entropy of an
    empty distribution is undefined, not zero.
    """

    policy = (
        molecule_type
        if isinstance(molecule_type, AlphabetPolicy)
        else policy_for(molecule_type or "protein")
        if molecule_type is not None
        else None
    )
    if policy is None:
        # Back-compatible behaviour: accept any non-gap, non-missing symbol.
        counts = Counter(v.upper() for v in values if v.upper() not in set("-.?NX"))
    else:
        counts = Counter(policy.normalize(v) for v in values if policy.is_residue(v))
    total = sum(counts.values())
    if total == 0:
        return None
    return -sum((c / total) * math.log2(c / total) for c in counts.values()) + 0.0


def entropy_from_counts(counts: np.ndarray) -> np.ndarray:
    """Row-wise Shannon entropy (bits) of a ``(L, K)`` count matrix; NaN where a row is empty."""

    totals = counts.sum(axis=1, keepdims=True).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(totals > 0, counts / totals, 0.0)
        logs = np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0.0)
        h = -(p * logs).sum(axis=1)
    h[totals[:, 0] == 0] = np.nan
    return h + 0.0
