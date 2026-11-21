"""Compositional heterogeneity and matched-pairs symmetry tests.

Most substitution models assume every sequence evolved under the same,
stationary composition. These diagnostics test that assumption on the alignment
alone (no tree):

* **RCV**, relative composition variability (Phillips & Penny 2003; as in PhyKIT):
  ``RCV = sum_c sum_i |x_ic - mean_c| / (N * L)``, with ``x_ic`` the count of
  residue ``c`` in sequence ``i``; per-sequence contributions (RCVT) sum to RCV.
* **Per-sequence χ² composition test** (as in IQ-TREE): each sequence's residue
  counts against the alignment-wide frequencies, ``df = k - 1``.
* **Jensen–Shannon divergence** of each sequence's composition from the alignment's.
* **Matched-pairs symmetry tests** (Bowker 1948; Stuart–Maxwell marginal symmetry;
  internal symmetry = Bowker − marginal), computed for each pair of sequences from
  their ``k x k`` divergence matrix, as in Naser-Khdour et al. (2019). A low p-value
  means the pair is unlikely to have evolved under the same stationary,
  reversible, homogeneous process.

Only canonical residues are counted.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from msatk.core.encoded import EncodedAlignment
from msatk.stats.partitions import _gammainc_upper
from msatk.stats.site_patterns import jensen_shannon

ALPHA = 0.05
DEFAULT_MAX_PAIRS = 20000


def chi2_sf(stat: float, df: float) -> float | None:
    if df <= 0 or not np.isfinite(stat):
        return None
    return float(_gammainc_upper(df / 2.0, max(stat, 0.0) / 2.0))


def composition_heterogeneity(
    enc: EncodedAlignment,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    counts = enc.sequence_counts.astype(np.float64)  # (N, k)
    n, k = counts.shape
    width = enc.width
    mean = counts.mean(axis=0)
    rcvt = np.abs(counts - mean).sum(axis=1) / (n * width) if n and width else np.zeros(n)
    totals = counts.sum(axis=0)
    freq = totals / totals.sum() if totals.sum() else np.zeros(k)
    present = freq > 0
    df = int(present.sum()) - 1
    rows = []
    jsd_values: list[float | None] = []
    failing = 0
    for i, seq_id in enumerate(enc.alignment.ids):
        m = counts[i].sum()
        if m > 0 and df > 0:
            expected = m * freq[present]
            stat = float(((counts[i, present] - expected) ** 2 / expected).sum())
            p = chi2_sf(stat, df)
            jsd = jensen_shannon(counts[i] / m, freq)
        else:
            stat, p, jsd = None, None, None  # type: ignore[assignment]
        jsd_values.append(jsd)
        fails = p is not None and p < ALPHA
        failing += int(fails)
        rows.append(
            {
                "sequence_id": seq_id,
                "rcvt": float(rcvt[i]),
                "composition_chi2": stat,
                "composition_df": df if stat is not None else None,
                "composition_p": p,
                "composition_test_failed": fails,
                "composition_jsd": jsd,
            }
        )
    tested = sum(r["composition_p"] is not None for r in rows)
    jsds = [j for j in jsd_values if j is not None]
    summary = {
        "rcv": float(rcvt.sum()),
        "max_rcvt": float(rcvt.max()) if n else None,
        "sequences_tested": tested,
        "sequences_failing_composition_test": failing,
        "fraction_failing_composition_test": failing / tested if tested else None,
        "alpha": ALPHA,
        "mean_composition_jsd": float(np.mean(jsds)) if jsds else None,
        "max_composition_jsd": float(np.max(jsds)) if jsds else None,
    }
    return summary, {"sequence_composition_tests": rows}


def _symmetry(m: np.ndarray) -> dict[str, Any]:
    """Bowker, Stuart–Maxwell marginal and internal symmetry for one k x k matrix."""

    k = m.shape[0]
    iu = np.triu_indices(k, 1)
    upper, lower = m[iu], m.T[iu]
    s = upper + lower
    use = s > 0
    bowker = float(((upper - lower)[use] ** 2 / s[use]).sum())
    df_b = int(use.sum())
    # Stuart–Maxwell on the first k-1 margins.
    d = (m.sum(axis=1) - m.sum(axis=0))[:-1]
    v = -(m + m.T)
    np.fill_diagonal(v, m.sum(axis=1) + m.sum(axis=0) - 2 * np.diag(m))
    v = v[:-1, :-1]
    df_m = int(np.linalg.matrix_rank(v)) if v.size else 0
    marginal = float(d @ np.linalg.pinv(v) @ d) if df_m else 0.0
    internal = bowker - marginal
    df_i = df_b - df_m
    return {
        "sites": int(m.sum()),
        "bowker": bowker,
        "bowker_df": df_b,
        "bowker_p": chi2_sf(bowker, df_b),
        "marginal": marginal,
        "marginal_df": df_m,
        "marginal_p": chi2_sf(marginal, df_m),
        "internal": internal,
        "internal_df": df_i,
        "internal_p": chi2_sf(internal, df_i) if df_i > 0 else None,
    }


def symmetry_tests(
    enc: EncodedAlignment, max_pairs: int = DEFAULT_MAX_PAIRS, seed: int = 0
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    n, k = enc.n, enc.k
    total = n * (n - 1) // 2
    if total <= max_pairs:
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        method = "exact"
    else:
        rng = np.random.default_rng(seed)
        chosen: set[tuple[int, int]] = set()
        while len(chosen) < max_pairs:
            i, j = sorted(rng.choice(n, 2, replace=False).tolist())
            chosen.add((i, j))
        pairs = sorted(chosen)
        method = "sampled"
    idx = enc.residue_index.astype(np.int64)
    rows = []
    for i, j in pairs:
        both = (idx[i] >= 0) & (idx[j] >= 0)
        m = np.bincount(idx[i][both] * k + idx[j][both], minlength=k * k).reshape(k, k)
        res = _symmetry(m.astype(np.float64))
        rows.append({"sequence_a": enc.alignment.ids[i], "sequence_b": enc.alignment.ids[j], **res})

    def frac(key: str) -> float | None:
        vals = [r[key] for r in rows if r[key] is not None]
        return float(np.mean([v < ALPHA for v in vals])) if vals else None

    def minimum(key: str) -> float | None:
        vals = [r[key] for r in rows if r[key] is not None]
        return float(min(vals)) if vals else None

    summary = {
        "pair_method": method,
        "pairs_tested": len(rows),
        "pairs_total": total,
        "sample_seed": seed if method == "sampled" else None,
        "alpha": ALPHA,
        "fraction_pairs_rejecting_bowker": frac("bowker_p"),
        "fraction_pairs_rejecting_marginal": frac("marginal_p"),
        "fraction_pairs_rejecting_internal": frac("internal_p"),
        "min_bowker_p": minimum("bowker_p"),
        "multiple_testing": "none applied; fractions are of uncorrected tests at alpha",
    }
    return summary, {"pair_symmetry_tests": rows}
