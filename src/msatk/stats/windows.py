"""Statistics in sliding windows along the alignment.

Global averages hide local problems: a well-aligned gene with one badly aligned
exon, or a hypervariable loop. Each window reports occupancy, gaps, variability
and entropy, and counts *window outlier sequences*: sequences whose identity to
the consensus inside the window has a robust z-score below ``-3.5`` (median/MAD
over sequences with at least half the window observed; the scaled MAD is floored at
0.02 so that a window where most sequences are identical still has a finite scale).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from msatk.core.alphabet import is_nucleotide
from msatk.core.encoded import STATE_GAP, EncodedAlignment
from msatk.stats.site import site_arrays

OUTLIER_Z = -3.5
MAD_FLOOR = 0.02  # identity units; avoids z = -inf when most rows are identical


def default_window(width: int, codon: bool) -> int:
    w = int(min(max(30, round(width / 20)), 300))
    if codon:
        w -= w % 3
    return max(min(w, width), 1)


def windowed_analysis(
    enc: EncodedAlignment, window: int | None = None, step: int | None = None
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    codon = enc.policy.molecule_type == "codon"
    width = enc.width
    window = window or default_window(width, codon)
    step = step or max(window // 2, 1)
    if codon and (window % 3 or step % 3):
        raise ValueError("Codon windows and steps must be multiples of 3.")
    arrays = site_arrays(enc)
    gap = (enc.state == STATE_GAP).mean(axis=0)
    occupancy = arrays["observed"] / max(enc.n, 1)
    entropy = arrays["entropy"]
    variable = arrays["variable"]
    informative = arrays["parsimony_informative"]
    counts = enc.site_counts
    nucleotide = is_nucleotide(enc.policy.molecule_type)
    if nucleotide:
        alphabet = enc.policy.canonical
        gc_counts = counts[:, alphabet.index("G")] + counts[:, alphabet.index("C")]
    cons = enc.consensus_index
    match = (enc.residue_index == cons[None, :]) & (cons[None, :] >= 0)
    obs = enc.residue_mask
    starts = list(range(0, max(width - window, 0) + 1, step))
    if starts and starts[-1] + window < width:
        starts.append(width - window if not codon else (width - window) - (width - window) % 3)
    rows = []
    for s in sorted(set(starts)):
        e = min(s + window, width)
        sl = slice(s, e)
        ent = entropy[sl]
        n_obs = obs[:, sl].sum(axis=1)
        ident = np.where(n_obs > 0, match[:, sl].sum(axis=1) / np.maximum(n_obs, 1), np.nan)
        enough = n_obs >= 0.5 * (e - s)
        vals = ident[enough]
        outliers = 0
        if vals.size >= 5:
            med = float(np.median(vals))
            mad = max(float(np.median(np.abs(vals - med))) * 1.4826, MAD_FLOOR)
            outliers = int(((vals - med) / mad < OUTLIER_Z).sum())
        row: dict[str, Any] = {
            "window_start": s + 1,
            "window_end": e,
            "sites": e - s,
            "occupancy": float(occupancy[sl].mean()),
            "gap_fraction": float(gap[sl].mean()),
            "variable_fraction": float(variable[sl].mean()),
            "parsimony_informative_fraction": float(informative[sl].mean()),
            "mean_entropy": float(np.nanmean(ent)) if np.isfinite(ent).any() else None,
            "mean_consensus_identity": float(np.nanmean(vals)) if vals.size else None,
            "outlier_sequences": outliers,
        }
        if nucleotide:
            tot = counts[sl].sum()
            row["gc_content"] = float(gc_counts[sl].sum() / tot) if tot else None
        rows.append(row)

    def top(key: str) -> dict[str, Any] | None:
        valid = [r for r in rows if r[key] is not None]
        if not valid:
            return None
        best = max(valid, key=lambda r: r[key])
        return {"start": best["window_start"], "end": best["window_end"], key: best[key]}

    summary = {
        "window": window,
        "step": step,
        "windows": len(rows),
        "most_variable_window": top("variable_fraction"),
        "gappiest_window": top("gap_fraction"),
        "highest_entropy_window": top("mean_entropy"),
        "windows_with_outlier_sequences": sum(r["outlier_sequences"] > 0 for r in rows),
        "outlier_z_threshold": OUTLIER_Z,
        "outlier_mad_floor": MAD_FLOOR,
    }
    return summary, {"windows": rows}
