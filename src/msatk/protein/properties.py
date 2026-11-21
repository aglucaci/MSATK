"""Physicochemical variation and substitution severity in protein alignments.

* Per column: standard deviation of Kyte–Doolittle hydropathy and residue volume
  (Zamyatnin 1972), variance of net charge, and the fraction of polar residues,
  over canonical residues. A column can vary in identity yet keep its character
  (low hydropathy SD) or change it.
* Substitutions relative to the column consensus, scored by BLOSUM62:
  *conservative* when the score is positive, *radical* when negative, *neutral*
  at 0 (cf. BMGE's substitution-aware entropy, Criscuolo & Gribaldo 2010).
* Residue-class changes against the consensus (hydrophobic, polar, positive,
  negative, special: G/P/C), counted per column.
* Per sequence: mean BLOSUM62 score against the consensus, and its robust z-score,
  flagging sequences with unusually radical replacements.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np

from msatk.constants import HYDROPATHY_KD
from msatk.core.encoded import EncodedAlignment
from msatk.encoding.tensors import _CHARGE, _POLAR, _VOLUME
from msatk.stats.substitution import _canonical_matrix

CLASSES = {
    **dict.fromkeys("AVILMFWY", "hydrophobic"),
    **dict.fromkeys("STNQH", "polar"),
    **dict.fromkeys("KR", "positive"),
    **dict.fromkeys("DE", "negative"),
    **dict.fromkeys("GPC", "special"),
}


def protein_property_analysis(
    enc: EncodedAlignment,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    canon = enc.policy.canonical
    hyd = np.array([HYDROPATHY_KD.get(a, np.nan) for a in canon])
    vol = np.array([_VOLUME.get(a, np.nan) for a in canon])
    chg = np.array([_CHARGE.get(a, 0.0) for a in canon])
    pol = np.array([a in _POLAR for a in canon], dtype=float)
    cls = np.array([CLASSES.get(a, "other") for a in canon])
    blosum = _canonical_matrix(canon)
    idx = enc.residue_index.astype(np.int64)
    obs = idx >= 0
    safe = np.where(obs, idx, 0)
    cons = enc.consensus_index.astype(np.int64)

    def col_stat(values: np.ndarray, fn: str) -> np.ndarray:
        v = np.where(obs, values[safe], np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN columns -> NaN
            if fn == "std":
                return np.nanstd(v, axis=0)
            if fn == "var":
                return np.nanvar(v, axis=0)
            return np.nanmean(v, axis=0)

    hyd_sd = col_stat(hyd, "std")
    vol_sd = col_stat(vol, "std")
    chg_var = col_stat(chg, "var")
    pol_frac = col_stat(pol, "mean")
    cons_ok = cons >= 0
    differs = obs & cons_ok[None, :] & (idx != cons[None, :])
    score = np.where(differs, blosum[safe, np.where(cons_ok, cons, 0)[None, :]], np.nan)
    conservative = differs & (score > 0)
    radical = differs & (score < 0)
    neutral = differs & (score == 0)
    class_change = differs & (cls[safe] != cls[np.where(cons_ok, cons, 0)][None, :])
    rows = []
    for j in range(enc.width):
        n_obs = int(obs[:, j].sum())
        n_sub = int(differs[:, j].sum())
        rows.append(
            {
                "site_index": j + 1,
                "hydropathy_sd": None if n_obs < 2 else float(hyd_sd[j]),
                "volume_sd": None if n_obs < 2 else float(vol_sd[j]),
                "charge_variance": None if n_obs < 2 else float(chg_var[j]),
                "polar_fraction": None if n_obs == 0 else float(pol_frac[j]),
                "substitutions_vs_consensus": n_sub,
                "conservative_substitutions": int(conservative[:, j].sum()),
                "radical_substitutions": int(radical[:, j].sum()),
                "class_changes": int(class_change[:, j].sum()),
                "radical_fraction": float(radical[:, j].sum() / n_sub) if n_sub else None,
            }
        )
    # Per sequence: mean BLOSUM62 to consensus over observed, consensus-defined cells.
    both = obs & cons_ok[None, :]
    seq_score = np.where(both, blosum[safe, np.where(cons_ok, cons, 0)[None, :]], np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean_score = np.nanmean(seq_score, axis=1)
    finite = np.isfinite(mean_score)
    med = float(np.median(mean_score[finite])) if finite.any() else np.nan
    mad = float(np.median(np.abs(mean_score[finite] - med))) * 1.4826 if finite.any() else 0.0
    seq_rows: list[dict[str, Any]] = []
    for i, seq_id in enumerate(enc.alignment.ids):
        m = mean_score[i]
        seq_rows.append(
            {
                "sequence_id": seq_id,
                "mean_blosum62_to_consensus": float(m) if np.isfinite(m) else None,
                "substitution_robust_z": float((m - med) / mad)
                if np.isfinite(m) and mad > 0
                else None,
                "radical_substitutions": int(radical[i].sum()),
                "conservative_substitutions": int(conservative[i].sum()),
            }
        )
    robust_z: list[float | None] = [r["substitution_robust_z"] for r in seq_rows]
    total_sub = int(differs.sum())
    summary = {
        "substitutions_vs_consensus": total_sub,
        "conservative_fraction": float(conservative.sum() / total_sub) if total_sub else None,
        "radical_fraction": float(radical.sum() / total_sub) if total_sub else None,
        "neutral_fraction": float(neutral.sum() / total_sub) if total_sub else None,
        "class_change_fraction": float(class_change.sum() / total_sub) if total_sub else None,
        "median_hydropathy_sd": float(np.nanmedian(hyd_sd)) if np.isfinite(hyd_sd).any() else None,
        "sequences_with_substitution_z_below_minus_3_5": sum(
            1 for z in robust_z if z is not None and z < -3.5
        ),
    }
    return summary, {"site_properties": rows, "sequence_substitution_scores": seq_rows}
