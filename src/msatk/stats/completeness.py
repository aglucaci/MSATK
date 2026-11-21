"""Alignment completeness (AliStat-style) and gap structure.

Completeness follows Wong et al. (2020, AliStat): a cell is *observed* when it
holds a canonical residue; gaps, missing data (``?``, ``N``/``X``), ambiguity codes
and padding are unobserved. Two alignments can have the same overall missing
fraction but very different pairwise overlap, which is what decides whether
pairwise distances and tree inference are reliable, so MSATK reports overlap
between every pair of sequences as well as the global figures.

Gap structure distinguishes isolated gaps from long runs, internal gaps from
terminal ones (unsequenced ends), and reports the longest contiguous regions
without gaps or without any unobserved cell.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from msatk.core.encoded import STATE_GAP, EncodedAlignment

PAIR_BLOCK = 512
DEFAULT_MIN_PAIR_COMPLETENESS = 0.5


def _q(values: np.ndarray, qs: tuple[float, ...] = (0.05, 0.25, 0.5, 0.75, 0.95)) -> dict[str, Any]:
    if values.size == 0:
        return {f"q{int(q * 100):02d}": None for q in qs}
    return {f"q{int(q * 100):02d}": float(np.quantile(values, q)) for q in qs}


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """``(start, length)`` of each True run in a 1-D boolean array."""

    if not mask.any():
        return []
    padded = np.concatenate([[False], mask, [False]]).astype(np.int8)
    diff = np.diff(padded)
    starts = np.flatnonzero(diff == 1)
    ends = np.flatnonzero(diff == -1)
    return list(zip(starts.tolist(), (ends - starts).tolist()))


def pairwise_completeness(enc: EncodedAlignment) -> np.ndarray:
    """``(N, N)`` fraction of all L columns observed in both sequences (AliStat C_ij)."""

    obs = enc.residue_mask.astype(np.float32)
    n, width = obs.shape
    out = np.zeros((n, n), dtype=np.float64)
    for start in range(0, n, PAIR_BLOCK):
        stop = min(n, start + PAIR_BLOCK)
        out[start:stop] = obs[start:stop] @ obs.T
    return np.rint(out) / max(width, 1)


def completeness_analysis(
    enc: EncodedAlignment,
    exact_pairs: bool = True,
    min_pair_completeness: float = DEFAULT_MIN_PAIR_COMPLETENESS,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]], np.ndarray | None]:
    """Summary, tables and (when ``exact_pairs``) the C_ij matrix."""

    obs = enc.residue_mask
    n, width = obs.shape
    row_c = obs.mean(axis=1) if width else np.zeros(n)
    col_c = obs.mean(axis=0) if n else np.zeros(width)
    summary: dict[str, Any] = {
        "alignment_completeness": float(obs.mean()) if obs.size else 0.0,
        "min_sequence_completeness": float(row_c.min()) if n else None,
        "min_site_completeness": float(col_c.min()) if width else None,
        "sequence_completeness_quantiles": _q(row_c),
        "site_completeness_quantiles": _q(col_c),
        "min_pair_completeness_threshold": min_pair_completeness,
    }
    seq_rows: list[dict[str, Any]] = []
    matrix: np.ndarray | None = None
    worst_partner: list[tuple[str | None, float | None]] = [(None, None)] * n
    if exact_pairs and n >= 2:
        matrix = pairwise_completeness(enc)
        iu = np.triu_indices(n, 1)
        pairs = matrix[iu]
        summary.update(
            {
                "pair_completeness_method": "exact",
                "mean_pair_completeness": float(pairs.mean()),
                "min_pair_completeness": float(pairs.min()),
                "pair_completeness_quantiles": _q(pairs),
                "pairs_below_threshold": int((pairs < min_pair_completeness).sum()),
                "pairs_without_overlap": int((pairs == 0).sum()),
                "fraction_pairs_below_threshold": float((pairs < min_pair_completeness).mean()),
                "min_shared_sites": int(round(float(pairs.min()) * width)),
            }
        )
        off = matrix.copy()
        np.fill_diagonal(off, np.inf)
        partner = off.argmin(axis=1)
        worst_partner = [
            (enc.alignment.ids[int(j)], float(off[i, j])) for i, j in enumerate(partner)
        ]
    else:
        summary["pair_completeness_method"] = (
            "not computed (fewer than 2 sequences)"
            if n < 2
            else "not computed (pairwise=none or above the exact-comparison limit)"
        )
    for i, seq_id in enumerate(enc.alignment.ids):
        partner_id, partner_c = worst_partner[i]
        seq_rows.append(
            {
                "sequence_id": seq_id,
                "sequence_completeness": float(row_c[i]),
                "least_overlapping_partner": partner_id,
                "min_pair_completeness": partner_c,
                "pairs_below_threshold": int(
                    (matrix[i] < min_pair_completeness).sum()
                    - (matrix[i, i] < min_pair_completeness)
                )
                if matrix is not None
                else None,
            }
        )
    return summary, {"sequence_completeness": seq_rows}, matrix


def gap_structure_analysis(
    enc: EncodedAlignment,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    """Gap runs, terminal versus internal gaps, contiguous blocks and gap-pattern diversity."""

    gap = enc.state == STATE_GAP
    obs = enc.residue_mask
    n, width = gap.shape
    run_lengths: list[int] = []
    seq_rows: list[dict[str, Any]] = []
    terminal_total = internal_total = 0
    fragmentation: list[float] = []
    for i, seq_id in enumerate(enc.alignment.ids):
        # Terminal gaps: gap cells before the first / after the last non-gap cell.
        non_gap = np.flatnonzero(~gap[i])
        if non_gap.size:
            first, last = int(non_gap[0]), int(non_gap[-1])
            terminal = int(gap[i, :first].sum() + gap[i, last + 1 :].sum())
            internal_mask = gap[i].copy()
            internal_mask[:first] = False
            internal_mask[last + 1 :] = False
        else:
            terminal, internal_mask = int(gap[i].sum()), np.zeros(width, dtype=bool)
        runs = _runs(internal_mask)
        lengths = [length for _, length in runs]
        run_lengths += lengths
        internal = int(internal_mask.sum())
        terminal_total += terminal
        internal_total += internal
        segments = _runs(obs[i])
        observed = int(obs[i].sum())
        longest_segment = max((length for _, length in segments), default=0)
        frag = 1.0 - longest_segment / observed if observed else None
        if frag is not None:
            fragmentation.append(frag)
        seq_rows.append(
            {
                "sequence_id": seq_id,
                "internal_gap_openings": len(runs),
                "internal_gap_extensions": internal - len(runs),
                "internal_gap_cells": internal,
                "terminal_gap_cells": terminal,
                "longest_internal_gap": max(lengths, default=0),
                "observed_segments": len(segments),
                "longest_observed_segment": longest_segment,
                "fragmentation": frag,
            }
        )
    # Gap-pattern outliers: Jaccard distance between gap-position sets.
    gap_f = gap.astype(np.float32)
    inter = gap_f @ gap_f.T if n <= 4000 else None
    if inter is not None and n >= 3:
        sizes = gap_f.sum(axis=1)
        union = sizes[:, None] + sizes[None, :] - inter
        with np.errstate(divide="ignore", invalid="ignore"):
            jacc = np.where(union > 0, inter / union, 1.0)
        np.fill_diagonal(jacc, np.nan)
        mean_sim = np.nanmean(jacc, axis=1)
        dist = 1.0 - mean_sim
        med = float(np.median(dist))
        mad = float(np.median(np.abs(dist - med))) * 1.4826
        iu = np.triu_indices(n, 1)
        concordance = float(np.nanmean(jacc[iu]))
        for i, row in enumerate(seq_rows):
            row["gap_pattern_distance"] = float(dist[i])
            row["gap_pattern_robust_z"] = float((dist[i] - med) / mad) if mad > 0 else None
    else:
        concordance = None
        for row in seq_rows:
            row["gap_pattern_distance"] = None
            row["gap_pattern_robust_z"] = None
    gap_free = np.asarray(~gap.any(axis=0))
    complete = np.asarray(obs.all(axis=0))
    free_runs = _runs(gap_free)
    complete_runs = _runs(complete)
    lengths_arr = np.array(run_lengths, dtype=float)
    hist: dict[int, int] = {}
    for length in run_lengths:
        hist[length] = hist.get(length, 0) + 1
    col_patterns = {gap[:, j].tobytes() for j in range(width) if gap[:, j].any()}
    row_patterns = {gap[i].tobytes() for i in range(n)}
    total_gaps = int(gap.sum())
    summary: dict[str, Any] = {
        "gap_openings": len(run_lengths),
        "gap_extensions": int(lengths_arr.sum() - len(run_lengths)) if run_lengths else 0,
        "internal_gap_cells": internal_total,
        "terminal_gap_cells": terminal_total,
        "terminal_gap_fraction": float(terminal_total / total_gaps) if total_gaps else None,
        "mean_gap_run_length": float(lengths_arr.mean()) if run_lengths else None,
        "median_gap_run_length": float(np.median(lengths_arr)) if run_lengths else None,
        "max_gap_run_length": int(lengths_arr.max()) if run_lengths else 0,
        "single_site_gap_runs": hist.get(1, 0),
        "longest_gap_free_block": max((length for _, length in free_runs), default=0),
        "gap_free_blocks": len(free_runs),
        "longest_complete_block": max((length for _, length in complete_runs), default=0),
        "complete_blocks": len(complete_runs),
        "distinct_column_gap_patterns": len(col_patterns),
        "distinct_sequence_gap_patterns": len(row_patterns),
        "mean_pairwise_gap_concordance": concordance,
        "fragmentation_index": float(np.mean(fragmentation)) if fragmentation else None,
    }
    tables = {
        "gap_structure_by_sequence": seq_rows,
        "gap_run_length_distribution": [
            {"run_length": k, "runs": v} for k, v in sorted(hist.items())
        ],
        "contiguous_blocks": [
            {"block_type": "gap_free", "start": s + 1, "end": s + length, "length": length}
            for s, length in free_runs
        ]
        + [
            {"block_type": "complete", "start": s + 1, "end": s + length, "length": length}
            for s, length in complete_runs
        ],
    }
    return summary, tables
