"""Per-sequence statistics."""

from __future__ import annotations

from typing import Any

import numpy as np

from msatk.core.alphabet import is_nucleotide
from msatk.core.encoded import STATE_AMBIGUOUS, STATE_MISSING, EncodedAlignment
from msatk.stats.pairwise import PairwiseResult


def _gap_runs(row: np.ndarray) -> tuple[int, int, int, int, int]:
    """Return (leading, trailing, internal_gap_cells, run_count, max_run) for a gap mask row."""

    width = row.size
    if width == 0:
        return 0, 0, 0, 0, 0
    if row.all():
        return width, 0, 0, 1, width
    first = int(np.argmin(row))
    last = width - 1 - int(np.argmin(row[::-1]))
    leading, trailing = first, width - 1 - last
    padded = np.concatenate(([False], row, [False])).astype(np.int8)
    diff = np.diff(padded)
    starts, ends = np.where(diff == 1)[0], np.where(diff == -1)[0]
    runs = ends - starts
    internal = int(row[first : last + 1].sum())
    return leading, trailing, internal, int(runs.size), int(runs.max()) if runs.size else 0


def _skew(a: float, b: float) -> float | None:
    return (a - b) / (a + b) if (a + b) else None


def per_sequence_stats_encoded(
    enc: EncodedAlignment, pairwise: PairwiseResult | None = None
) -> list[dict[str, Any]]:
    width = enc.width
    alphabet = enc.policy.canonical
    nucleotide = is_nucleotide(enc.policy.molecule_type)
    residues = enc.residue_mask
    gaps = enc.gap_mask
    missing = (enc.state == STATE_MISSING) & ~enc.pad_mask
    ambiguous = enc.state == STATE_AMBIGUOUS
    consensus = enc.consensus_index
    consensus_defined = consensus >= 0
    seq_counts = enc.sequence_counts

    # Composition deviation: total-variation distance from the pooled residue composition.
    pooled = seq_counts.sum(axis=0).astype(float)
    pooled = pooled / pooled.sum() if pooled.sum() else pooled

    strings = enc.row_strings()
    first_seen: dict[str, str] = {}

    identity = pairwise.identity if pairwise is not None and pairwise.is_exact else None
    rows: list[dict[str, Any]] = []
    for i, record in enumerate(enc.alignment.records):
        leading, trailing, internal, run_count, max_run = _gap_runs(gaps[i])
        residue_count = int(residues[i].sum())
        comparable = residues[i] & consensus_defined
        n_comp = int(comparable.sum())
        agree = int((enc.residue_index[i][comparable] == consensus[comparable]).sum())
        cons_identity = agree / n_comp if n_comp else None
        counts = seq_counts[i].astype(float)
        comp_dist = (
            float(0.5 * np.abs(counts / counts.sum() - pooled).sum()) if counts.sum() else None
        )
        non_gap = np.where(~gaps[i] & ~enc.pad_mask[i])[0]
        key = strings[i]
        duplicate_of = first_seen.get(key)
        first_seen.setdefault(key, record.id)
        row: dict[str, Any] = {
            "row_index": i + 1,
            "sequence_id": record.id,
            "original_length": int(enc.original_lengths[i]),
            "aligned_width": width,
            "padded_cells": int(enc.pad_mask[i].sum()),
            "ungapped_length": int((~gaps[i] & ~enc.pad_mask[i]).sum()),
            "residue_count": residue_count,
            "span_start": int(non_gap[0] + 1) if non_gap.size else None,
            "span_end": int(non_gap[-1] + 1) if non_gap.size else None,
            "span_length": int(non_gap[-1] - non_gap[0] + 1) if non_gap.size else 0,
            "gap_count": int(gaps[i].sum()),
            "gap_fraction": float(gaps[i].sum() / width) if width else 0.0,
            "leading_gaps": leading,
            "trailing_gaps": trailing,
            "internal_gap_count": internal,
            "gap_run_count": run_count,
            "max_gap_run": max_run,
            "missing_count": int(missing[i].sum()),
            "ambiguous_count": int(ambiguous[i].sum()),
            "ambiguous_fraction": float(ambiguous[i].sum() / width) if width else 0.0,
            "occupancy": residue_count / width if width else 0.0,
            "unique_residue_count": int((seq_counts[i] > 0).sum()),
            "consensus_comparable_sites": n_comp,
            "consensus_identity": cons_identity,
            "composition_distance": comp_dist,
            "duplicate_of": duplicate_of,
        }
        if nucleotide:
            c = dict(zip(alphabet, seq_counts[i].tolist()))
            t = c.get("T", c.get("U", 0))
            gc_total = c["G"] + c["C"]
            row["gc_content"] = gc_total / residue_count if residue_count else None
            row["gc_skew"] = _skew(c["G"], c["C"])
            row["at_skew"] = _skew(c["A"], t)
        if identity is not None:
            others = np.delete(identity[i], i)
            defined = others[~np.isnan(others)]
            row["mean_identity_to_others"] = float(defined.mean()) if defined.size else None
            if defined.size:
                masked = np.where(np.isnan(identity[i]), -1.0, identity[i])
                masked[i] = -1.0
                j = int(masked.argmax())
                row["nearest_neighbor"] = enc.alignment.ids[j]
                row["nearest_neighbor_identity"] = float(identity[i, j])
            else:
                row["nearest_neighbor"] = None
                row["nearest_neighbor_identity"] = None
        else:
            row["mean_identity_to_others"] = None
            row["nearest_neighbor"] = None
            row["nearest_neighbor_identity"] = None
        row["outlier_score"] = None if cons_identity is None else 1.0 - cons_identity
        rows.append(row)
    _add_robust_z(rows, "consensus_identity", "consensus_identity_robust_z")
    _add_robust_z(rows, "composition_distance", "composition_robust_z")
    return rows


def _add_robust_z(rows: list[dict[str, Any]], key: str, out: str) -> None:
    """Median/MAD-based z-score (scaled by 1.4826); ``None`` when undefined."""

    values = np.array([np.nan if r[key] is None else r[key] for r in rows], dtype=float)
    finite = values[~np.isnan(values)]
    if finite.size < 3:
        for r in rows:
            r[out] = None
        return
    med = np.median(finite)
    scale = 1.4826 * np.median(np.abs(finite - med))
    if scale == 0:
        # Fall back to the mean absolute deviation (consistent for the normal distribution).
        scale = 1.2533 * np.mean(np.abs(finite - med))
    for r, v in zip(rows, values):
        if np.isnan(v):
            r[out] = None
        else:
            r[out] = 0.0 if scale == 0 else float((v - med) / scale)
