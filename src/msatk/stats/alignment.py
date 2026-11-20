"""Alignment-level summary statistics and the combined statistics bundle."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from statistics import mean, median, pstdev
from typing import Any

import numpy as np

from msatk.core.alphabet import is_nucleotide
from msatk.core.detect import alignment_detection_summary
from msatk.core.encoded import (
    STATE_AMBIGUOUS,
    STATE_GAP,
    STATE_INVALID,
    STATE_MISSING,
    STATE_RESIDUE,
    STATE_STOP,
    EncodedAlignment,
    encode_alignment,
)
from msatk.models import Alignment
from msatk.stats.diversity import diversity_summary
from msatk.stats.pairwise import PairwiseResult, compute_pairwise
from msatk.stats.sequence import per_sequence_stats_encoded
from msatk.stats.site import per_site_stats_encoded


def _mean(values: list[float]) -> float | None:
    return mean(values) if values else None


def _median(values: list[float]) -> float | None:
    return median(values) if values else None


def composition_rows(enc: EncodedAlignment) -> list[dict[str, Any]]:
    """Counts of every symbol, with its state, over the whole alignment (padding excluded)."""

    from msatk.core.encoded import STATE_NAMES

    real = ~enc.pad_mask
    keys = enc.chars[real].astype(np.int32) * 8 + enc.state[real].astype(np.int32)
    uniq, cnt = np.unique(keys, return_counts=True)
    counter = {(int(k) // 8, int(k) % 8): int(c) for k, c in zip(uniq, cnt)}
    residue_total = int(enc.residue_mask.sum())
    total = int(real.sum())
    rows = []
    for (code, state), count in sorted(counter.items()):
        rows.append(
            {
                "character": chr(code),
                "state": STATE_NAMES[state],
                "count": int(count),
                "fraction": count / total if total else 0.0,
                "fraction_of_residues": count / residue_total
                if state == STATE_RESIDUE and residue_total
                else None,
            }
        )
    return rows


def gap_block_rows(enc: EncodedAlignment, min_length: int = 1) -> list[dict[str, Any]]:
    """Contiguous gap runs per sequence (indel blocks), with terminal/internal classification."""

    rows: list[dict[str, Any]] = []
    for i, record in enumerate(enc.alignment.records):
        mask = enc.gap_mask[i]
        padded = np.concatenate(([False], mask, [False])).astype(np.int8)
        diff = np.diff(padded)
        for s, e in zip(np.where(diff == 1)[0], np.where(diff == -1)[0]):
            if e - s < min_length:
                continue
            rows.append(
                {
                    "sequence_id": record.id,
                    "start": int(s + 1),
                    "end": int(e),
                    "length": int(e - s),
                    "position": "leading"
                    if s == 0
                    else ("trailing" if e == enc.width else "internal"),
                }
            )
    return rows


def alignment_summary_encoded(
    enc: EncodedAlignment,
    site_rows: list[dict[str, Any]],
    seq_rows: list[dict[str, Any]],
    pairwise: PairwiseResult | None,
) -> dict[str, Any]:
    alignment = enc.alignment
    molecule_type = enc.policy.molecule_type
    n, width = enc.n, enc.width
    cells = n * width
    real_cells = cells - alignment.padded_cells
    state = enc.state
    lengths = [row["ungapped_length"] for row in seq_rows]
    id_counts = Counter(alignment.ids)
    observed_sites = [row for row in site_rows if row["site_class"] != "unobserved"]
    entropies = [row["entropy"] for row in observed_sites if row["entropy"] is not None]
    norm_entropies = [
        row["normalized_entropy"] for row in observed_sites if row["normalized_entropy"] is not None
    ]
    cons_ids = [
        row["consensus_identity"] for row in seq_rows if row["consensus_identity"] is not None
    ]
    gaps = int((state == STATE_GAP).sum())
    summary: dict[str, Any] = {
        "tool": "MSATK",
        "number_of_sequences": n,
        "alignment_length": width,
        "total_cells": cells,
        "rectangular_alignment": alignment.is_rectangular,
        "sequence_lengths_consistent": alignment.is_rectangular,
        "padded_cells": alignment.padded_cells,
        "detected_format": alignment.fmt,
        "molecule_type": molecule_type,
        **alignment_detection_summary(alignment, molecule_type),
        "duplicate_id_count": sum(1 for c in id_counts.values() if c > 1),
        "unique_sequence_count": len(
            {row["duplicate_of"] or row["sequence_id"] for row in seq_rows}
        ),
        "duplicate_sequence_count": sum(1 for row in seq_rows if row["duplicate_of"]),
        "residue_cells": int(enc.residue_mask.sum()),
        "ambiguous_cells": int((state == STATE_AMBIGUOUS).sum()),
        "missing_cells": int(((state == STATE_MISSING) & ~enc.pad_mask).sum()),
        "stop_cells": int((state == STATE_STOP).sum()),
        "invalid_cells": int((state == STATE_INVALID).sum()),
        "total_gaps": gaps,
        "gap_fraction": gaps / cells if cells else 0.0,
        "occupancy": float(enc.residue_mask.sum() / cells) if cells else 0.0,
        "ambiguous_fraction": float((state == STATE_AMBIGUOUS).sum() / real_cells)
        if real_cells
        else 0.0,
        "mean_sequence_length": _mean(lengths),
        "median_sequence_length": _median(lengths),
        "min_sequence_length": min(lengths) if lengths else None,
        "max_sequence_length": max(lengths) if lengths else None,
        "sequence_length_sd": pstdev(lengths) if len(lengths) > 1 else 0.0,
        "constant_sites": sum(row["site_class"] == "constant" for row in site_rows),
        "conserved_sites": sum(row["site_class"] == "constant" for row in site_rows),
        "variable_sites": sum(row["site_class"] == "variable" for row in site_rows),
        "unobserved_sites": sum(row["site_class"] == "unobserved" for row in site_rows),
        "parsimony_informative_sites": sum(bool(row["parsimony_informative"]) for row in site_rows),
        "singleton_sites": sum(bool(row["singleton"]) for row in site_rows),
        "gap_free_sites": sum(row["gap_count"] == 0 for row in site_rows),
        "complete_sites": int(enc.residue_mask.all(axis=0).sum()),
        "entropy_mean": _mean(entropies),
        "entropy_median": _median(entropies),
        "normalized_entropy_mean": _mean(norm_entropies),
        "mean_consensus_identity": _mean(cons_ids),
    }
    if pairwise is not None:
        summary.update(pairwise.summary())
    else:
        summary.update(
            {
                "pairwise_method": "none",
                "mean_pairwise_identity": None,
                "mean_pairwise_distance": None,
            }
        )
    if is_nucleotide(molecule_type):
        counts = enc.site_counts.sum(axis=0)
        alphabet = enc.policy.canonical
        c = dict(zip(alphabet, counts.tolist()))
        t = c.get("T", c.get("U", 0))
        total = int(counts.sum())
        summary["gc_content"] = (c["G"] + c["C"]) / total if total else None
        summary["gc_skew"] = (c["G"] - c["C"]) / (c["G"] + c["C"]) if c["G"] + c["C"] else None
        summary["at_skew"] = (c["A"] - t) / (c["A"] + t) if c["A"] + t else None
        summary.update(diversity_summary(enc))
    else:
        summary["gc_content"] = None
    return summary


@dataclass
class StatisticsBundle:
    """All statistics computed for an alignment at once (sharing one encoding)."""

    encoded: EncodedAlignment
    summary: dict[str, Any]
    sites: list[dict[str, Any]]
    sequences: list[dict[str, Any]]
    pairwise: PairwiseResult | None
    composition: list[dict[str, Any]] = field(default_factory=list)


def compute_statistics(
    alignment: Alignment,
    molecule_type: str,
    pairwise: str = "auto",
    max_exact_sequences: int = 1000,
    sample_pairs: int = 20000,
    seed: int = 0,
) -> StatisticsBundle:
    enc = encode_alignment(alignment, molecule_type)
    pw = compute_pairwise(
        enc,
        method=pairwise,
        max_exact_sequences=max_exact_sequences,
        sample_pairs=sample_pairs,
        seed=seed,
    )
    sites = per_site_stats_encoded(enc)
    sequences = per_sequence_stats_encoded(enc, pw)
    summary = alignment_summary_encoded(enc, sites, sequences, pw)
    return StatisticsBundle(enc, summary, sites, sequences, pw, composition_rows(enc))
