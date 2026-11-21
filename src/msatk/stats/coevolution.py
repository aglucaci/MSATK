"""Mutual information between alignment columns, with APC correction.

For columns X and Y,

    I(X;Y) = sum_{x,y} p(x,y) log2( p(x,y) / (p(x) p(y)) ).

Raw MI is inflated by shared ancestry, redundant sequences, small samples and
column entropy. MSATK therefore applies, and records:

* **sequence weights** ``w_i = 1/n_i`` at 80 % identity (see
  :mod:`msatk.stats.redundancy`) when exact pairwise identities are available;
* a **pseudocount** ``lambda`` mixing observed frequencies with a uniform
  distribution: ``p = (1 - lambda) f + lambda / q`` (joint: ``lambda / q^2``);
* the **average product correction** (Dunn et al. 2008):
  ``APC(i,j) = MI(i,j) - MI(i,.) MI(.,j) / MI(.,.)``;
* normalised MI ``2 I / (H_X + H_Y)`` (symmetric uncertainty, 0..1).

Each pair uses only sequences with a canonical residue in both columns. To stay
tractable the analysis runs on candidate columns: variable, at least
``min_occupancy`` observed, the ``max_sites`` most entropic. This is a
screening statistic, not a contact predictor: direct-coupling models
(plmDCA, GREMLIN) are out of scope.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from msatk.core.encoded import EncodedAlignment

DEFAULT_MAX_SITES = 200
DEFAULT_MIN_OCCUPANCY = 0.5
DEFAULT_PSEUDOCOUNT = 0.5
DEFAULT_TOP_PAIRS = 200
DEFAULT_MIN_SEPARATION = 1


def coevolution_analysis(
    enc: EncodedAlignment,
    weights: np.ndarray | None = None,
    max_sites: int = DEFAULT_MAX_SITES,
    min_occupancy: float = DEFAULT_MIN_OCCUPANCY,
    pseudocount: float = DEFAULT_PSEUDOCOUNT,
    top_pairs: int = DEFAULT_TOP_PAIRS,
    min_separation: int = DEFAULT_MIN_SEPARATION,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    n, k = enc.n, enc.k
    w = np.ones(n) if weights is None else np.asarray(weights, dtype=np.float64)
    counts = enc.site_counts
    observed = counts.sum(axis=1)
    occ = observed / max(n, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(observed[:, None] > 0, counts / np.maximum(observed[:, None], 1), 0.0)
        h = -(np.where(p > 0, p * np.log2(np.where(p > 0, p, 1)), 0)).sum(axis=1)
    candidates = np.flatnonzero(((counts > 0).sum(axis=1) >= 2) & (occ >= min_occupancy))
    if candidates.size > max_sites:
        order = np.argsort(-h[candidates], kind="stable")[:max_sites]
        candidates = np.sort(candidates[order])
    s = candidates.size
    base = {
        "candidate_sites": int(s),
        "max_sites": max_sites,
        "min_occupancy": min_occupancy,
        "pseudocount": pseudocount,
        "weighting": "1/n_i at 80% identity" if weights is not None else "none (uniform)",
        "effective_sequences": float(w.sum()),
        "min_separation": min_separation,
    }
    if s < 2:
        return base | {"pairs": 0, "note": "fewer than 2 candidate columns"}, {
            "coevolution_pairs": []
        }
    # One-hot (N, S*k) design; weighted joint counts computed in blocks of sites.
    idx = enc.residue_index[:, candidates].astype(np.int64)
    onehot = np.zeros((n, s, k), dtype=np.float64)
    rows, cols = np.nonzero(idx >= 0)
    onehot[rows, cols, idx[rows, cols]] = 1.0
    a = onehot.reshape(n, s * k)
    aw = a * w[:, None]
    lam = pseudocount
    mi = np.zeros((s, s))
    hi = np.zeros((s, s))
    hj = np.zeros((s, s))
    hij = np.zeros((s, s))
    total = np.zeros((s, s))
    block = max(1, 400_000 // max(s * k * k, 1))
    for b0 in range(0, s, block):
        b1 = min(s, b0 + block)
        joint = aw[:, b0 * k : b1 * k].T @ a  # (b*k, S*k) weighted joint counts
        joint = joint.reshape(b1 - b0, k, s, k).transpose(0, 2, 1, 3)  # (b, S, k, k)
        tot = joint.sum(axis=(2, 3))
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(tot[..., None, None] > 0, joint / tot[..., None, None], 0.0)
            pij = (1 - lam) * f + lam / (k * k)
            pi = pij.sum(axis=3)
            pj = pij.sum(axis=2)
            ratio = pij / (pi[..., :, None] * pj[..., None, :])
            mi[b0:b1] = np.where(pij > 0, pij * np.log2(ratio), 0.0).sum(axis=(2, 3))
            hi[b0:b1] = -np.where(pi > 0, pi * np.log2(pi), 0.0).sum(axis=2)
            hj[b0:b1] = -np.where(pj > 0, pj * np.log2(pj), 0.0).sum(axis=2)
            hij[b0:b1] = -np.where(pij > 0, pij * np.log2(pij), 0.0).sum(axis=(2, 3))
        total[b0:b1] = tot
    mi = np.where(total > 0, mi, np.nan)
    np.fill_diagonal(mi, np.nan)
    row_mean = np.nanmean(mi, axis=1)
    all_mean = float(np.nanmean(mi))
    apc = mi - np.outer(row_mean, row_mean) / all_mean if all_mean > 0 else mi
    with np.errstate(divide="ignore", invalid="ignore"):
        nmi = np.where(hi + hj > 0, 2 * mi / (hi + hj), np.nan)
    iu = np.triu_indices(s, 1)
    sep = np.abs(candidates[iu[1]] - candidates[iu[0]])
    keep = (sep >= min_separation) & np.isfinite(mi[iu])
    order = np.argsort(-apc[iu][keep], kind="stable")[:top_pairs]
    ii, jj = iu[0][keep][order], iu[1][keep][order]
    out_rows = [
        {
            "rank": r + 1,
            "site_i": int(candidates[i]) + 1,
            "site_j": int(candidates[j]) + 1,
            "mi": float(mi[i, j]),
            "mi_apc": float(apc[i, j]),
            "normalized_mi": float(nmi[i, j]) if np.isfinite(nmi[i, j]) else None,
            "joint_entropy": float(hij[i, j]),
            "conditional_entropy_j_given_i": float(hij[i, j] - hi[i, j]),
            "effective_pair_weight": float(total[i, j]),
        }
        for r, (i, j) in enumerate(zip(ii.tolist(), jj.tolist()))
    ]
    vals = apc[iu][keep]
    summary = base | {
        "pairs": int(keep.sum()),
        "mean_mi": float(np.nanmean(mi[iu][keep])) if keep.any() else None,
        "max_mi_apc": float(np.nanmax(vals)) if vals.size else None,
        "top_pair": [out_rows[0]["site_i"], out_rows[0]["site_j"]] if out_rows else None,
        "pairs_reported": len(out_rows),
    }
    return summary, {"coevolution_pairs": out_rows}
