"""Per-partition statistics for multi-locus / partitioned alignments.

Each partition (gene, codon position, ...) is analysed as its own alignment with
the same engine, and compositional heterogeneity across sequences is tested with
Pearson's χ² on the sequence × residue count table (as in PAUP*/IQ-TREE's
composition test). Partitions keep their own boundaries; nothing is concatenated
silently.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from msatk.core.encoded import EncodedAlignment, encode_alignment
from msatk.io.readers import Partition
from msatk.models import Alignment, SequenceRecord


def _gammainc_upper(a: float, x: float) -> float:
    """Regularised upper incomplete gamma Q(a, x) (Numerical Recipes series / CF)."""

    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1:
        term = total = 1.0 / a
        ap = a
        for _ in range(1000):
            ap += 1
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-14:
                break
        return max(0.0, 1.0 - total * math.exp(-x + a * math.log(x) - gln))
    b = x + 1 - a
    c = 1 / 1e-300
    d = 1 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-14:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def composition_chi2(enc: EncodedAlignment) -> dict[str, Any]:
    """Pearson χ² test of compositional homogeneity across sequences."""

    counts = enc.sequence_counts.astype(np.float64)
    counts = counts[counts.sum(axis=1) > 0]
    counts = counts[:, counts.sum(axis=0) > 0]
    if counts.shape[0] < 2 or counts.shape[1] < 2:
        return {"composition_chi2": None, "composition_chi2_df": None, "composition_chi2_p": None}
    expected = counts.sum(axis=1, keepdims=True) * counts.sum(axis=0, keepdims=True) / counts.sum()
    chi2 = float(((counts - expected) ** 2 / expected).sum())
    df = (counts.shape[0] - 1) * (counts.shape[1] - 1)
    return {
        "composition_chi2": chi2,
        "composition_chi2_df": df,
        "composition_chi2_p": _gammainc_upper(df / 2, chi2 / 2),
    }


def partition_rows(
    alignment: Alignment, molecule_type: str, partitions: list[Partition]
) -> list[dict[str, Any]]:
    from msatk.stats.alignment import compute_statistics

    rows: list[dict[str, Any]] = []
    width = alignment.length
    for part in partitions:
        sites = [s for s in part.sites if 1 <= s <= width]
        row: dict[str, Any] = {
            "partition": part.name,
            "datatype": part.datatype or None,
            "sites_declared": len(part.sites),
            "sites_in_alignment": len(sites),
        }
        if sites:
            padded = alignment.padded_sequences()
            part_aln = Alignment(
                tuple(
                    SequenceRecord(r.id, "".join(seq[k - 1] for k in sites))
                    for r, seq in zip(alignment.records, padded)
                ),
                source=alignment.source,
                fmt=alignment.fmt,
            )
            ptype = "dna" if molecule_type == "codon" else molecule_type
            s = compute_statistics(part_aln, ptype, pairwise="none").summary
            enc = encode_alignment(part_aln, ptype)
            empty = sum(1 for r in part_aln.records if not any(c not in "-.?" for c in r.sequence))
            row.update(
                {
                    "occupancy": s["occupancy"],
                    "gap_fraction": s["gap_fraction"],
                    "sequences_without_data": empty,
                    "variable_sites": s["variable_sites"],
                    "parsimony_informative_sites": s["parsimony_informative_sites"],
                    "entropy_mean": s["entropy_mean"],
                    "gc_content": s.get("gc_content"),
                    **composition_chi2(enc),
                }
            )
        rows.append(row)
    return rows


def partition_summary(
    partitions: list[Partition], rows: list[dict[str, Any]], width: int
) -> dict[str, Any]:
    covered = len({s for part in partitions for s in part.sites if 1 <= s <= width})
    overlapping = sum(len(part.sites) for part in partitions) - len(
        {s for part in partitions for s in part.sites}
    )
    p: list[float] = [
        float(r["composition_chi2_p"]) for r in rows if r.get("composition_chi2_p") is not None
    ]
    return {
        "partitions": len(rows),
        "partition_sites_assigned": covered,
        "partition_sites_unassigned": max(width - covered, 0),
        "partition_sites_overlapping": overlapping,
        "partitions_failing_composition_test": sum(1 for v in p if float(v) < 0.05),
    }
