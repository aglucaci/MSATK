"""Pairwise sequence identity and p-distance with explicit cost policies.

Exact all-vs-all comparison is O(N^2 * L). MSATK computes it with blocked
matrix products (``sum_c X_c X_c^T``) over one-hot residue indicators, and
switches to a deterministic random sample of pairs above a size threshold so
that default profiling of large alignments stays tractable.

Identity between two sequences is ``matches / comparable_sites`` where a site is
comparable only when *both* sequences carry a canonical residue. When two
sequences share no comparable site their identity is **undefined** (NaN / null),
never zero.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from msatk.core.encoded import EncodedAlignment

PAIRWISE_METHODS = ("auto", "exact", "sampled", "none")
DEFAULT_MAX_EXACT_SEQUENCES = 1000
DEFAULT_SAMPLE_PAIRS = 20000


@dataclass
class PairwiseResult:
    """Pairwise comparison results.

    ``identity``/``comparable``/``mismatches`` are full ``(N, N)`` matrices for the
    ``exact`` method and ``None`` otherwise; ``pair_identities`` always holds the
    identities of every evaluated pair (``i < j``).
    """

    method: str
    ids: list[str]
    pair_identities: np.ndarray
    pairs_evaluated: int
    pairs_total: int
    pairs_without_comparable_sites: int
    identity: np.ndarray | None = None
    comparable: np.ndarray | None = None
    mismatches: np.ndarray | None = None
    seed: int | None = None

    @property
    def is_exact(self) -> bool:
        return self.identity is not None

    def distance(self) -> np.ndarray | None:
        return None if self.identity is None else 1.0 - self.identity

    def summary(self) -> dict[str, object]:
        values = self.pair_identities[~np.isnan(self.pair_identities)]
        defined = values.size > 0
        return {
            "pairwise_method": self.method,
            "pairwise_pairs_total": self.pairs_total,
            "pairwise_pairs_evaluated": self.pairs_evaluated,
            "pairwise_pairs_without_comparable_sites": self.pairs_without_comparable_sites,
            "pairwise_sample_seed": self.seed,
            "mean_pairwise_identity": float(values.mean()) if defined else None,
            "median_pairwise_identity": float(np.median(values)) if defined else None,
            "min_pairwise_identity": float(values.min()) if defined else None,
            "max_pairwise_identity": float(values.max()) if defined else None,
            "mean_pairwise_distance": float(1.0 - values.mean()) if defined else None,
        }

    def matrix_rows(self, which: str = "identity") -> list[list[float | None]]:
        matrix = self.identity if which == "identity" else self.distance()
        if matrix is None:
            raise ValueError("Full matrices are only available for exact pairwise comparison.")
        return [[None if np.isnan(v) else float(v) for v in row] for row in matrix]

    def long_table(self) -> list[dict[str, object]]:
        """One record per pair (exact mode only)."""

        if self.identity is None or self.comparable is None or self.mismatches is None:
            return []
        rows: list[dict[str, object]] = []
        n = len(self.ids)
        for i in range(n):
            for j in range(i + 1, n):
                ident = self.identity[i, j]
                rows.append(
                    {
                        "sequence_a": self.ids[i],
                        "sequence_b": self.ids[j],
                        "comparable_sites": int(self.comparable[i, j]),
                        "mismatches": int(self.mismatches[i, j]),
                        "identity": None if np.isnan(ident) else float(ident),
                        "p_distance": None if np.isnan(ident) else float(1.0 - ident),
                    }
                )
        return rows


def _indicators(enc: EncodedAlignment, rows: np.ndarray | slice) -> list[np.ndarray]:
    idx = enc.residue_index[rows]
    return [(idx == c).astype(np.float32) for c in range(enc.k)]


def exact_pairwise(enc: EncodedAlignment, block: int = 512) -> PairwiseResult:
    n = enc.n
    mask = enc.residue_mask.astype(np.float32)
    comparable = np.zeros((n, n), dtype=np.float64)
    matches = np.zeros((n, n), dtype=np.float64)
    full = _indicators(enc, slice(None))
    for start in range(0, n, block):
        stop = min(n, start + block)
        comparable[start:stop] = mask[start:stop] @ mask.T
        acc = np.zeros((stop - start, n), dtype=np.float64)
        for channel in full:
            acc += channel[start:stop] @ channel.T
        matches[start:stop] = acc
    comparable = np.rint(comparable)
    matches = np.rint(matches)
    with np.errstate(divide="ignore", invalid="ignore"):
        identity = np.where(comparable > 0, matches / comparable, np.nan)
    iu = np.triu_indices(n, k=1)
    pair_values = identity[iu]
    return PairwiseResult(
        method="exact",
        ids=enc.alignment.ids,
        pair_identities=pair_values,
        pairs_evaluated=int(pair_values.size),
        pairs_total=n * (n - 1) // 2,
        pairs_without_comparable_sites=int(np.isnan(pair_values).sum()),
        identity=identity,
        comparable=comparable.astype(np.int64),
        mismatches=(comparable - matches).astype(np.int64),
    )


def sampled_pairwise(
    enc: EncodedAlignment, n_pairs: int = DEFAULT_SAMPLE_PAIRS, seed: int = 0, block: int = 2048
) -> PairwiseResult:
    n = enc.n
    total = n * (n - 1) // 2
    rng = np.random.default_rng(seed)
    if total <= n_pairs:
        a, b = np.triu_indices(n, k=1)
    else:
        a = rng.integers(0, n, size=n_pairs * 2)
        b = rng.integers(0, n, size=n_pairs * 2)
        keep = a != b
        a, b = np.minimum(a[keep], b[keep])[:n_pairs], np.maximum(a[keep], b[keep])[:n_pairs]
    values = np.empty(a.size, dtype=float)
    idx = enc.residue_index
    mask = enc.residue_mask
    for start in range(0, a.size, block):
        sa, sb = a[start : start + block], b[start : start + block]
        both = mask[sa] & mask[sb]
        comp = both.sum(axis=1)
        same = ((idx[sa] == idx[sb]) & both).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            values[start : start + block] = np.where(comp > 0, same / comp, np.nan)
    return PairwiseResult(
        method="sampled" if total > n_pairs else "all_pairs",
        ids=enc.alignment.ids,
        pair_identities=values,
        pairs_evaluated=int(values.size),
        pairs_total=total,
        pairs_without_comparable_sites=int(np.isnan(values).sum()),
        seed=seed if total > n_pairs else None,
    )


def compute_pairwise(
    enc: EncodedAlignment,
    method: str = "auto",
    max_exact_sequences: int = DEFAULT_MAX_EXACT_SEQUENCES,
    sample_pairs: int = DEFAULT_SAMPLE_PAIRS,
    seed: int = 0,
) -> PairwiseResult | None:
    """Compute pairwise identity according to ``method``.

    ``auto`` → exact when ``N <= max_exact_sequences`` else sampled; ``none`` → ``None``.
    """

    if method not in PAIRWISE_METHODS:
        raise ValueError(f"Unknown pairwise method {method!r}; choose from {PAIRWISE_METHODS}.")
    if method == "none" or enc.n < 2:
        return None
    if method == "exact" or (method == "auto" and enc.n <= max_exact_sequences):
        return exact_pairwise(enc)
    return sampled_pairwise(enc, n_pairs=sample_pairs, seed=seed)


def pairwise_identity(seq_a: str, seq_b: str, molecule_type: str = "dna") -> float | None:
    """Identity of two aligned strings over mutually comparable canonical sites (``None`` if none)."""

    from msatk.core.alphabet import policy_for

    policy = policy_for(molecule_type)
    matches = comparable = 0
    for a, b in zip(seq_a, seq_b):
        if policy.is_residue(a) and policy.is_residue(b):
            comparable += 1
            matches += int(policy.normalize(a) == policy.normalize(b))
    return matches / comparable if comparable else None
