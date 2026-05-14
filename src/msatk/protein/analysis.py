"""Protein composition and residue-class summaries."""

from __future__ import annotations

from collections import Counter
from statistics import mean

from msatk.constants import AA_CANONICAL, HYDROPATHY_KD, MISSING_CHARS, RESIDUE_CLASSES
from msatk.models import Alignment


def amino_acid_composition(alignment: Alignment) -> list[dict[str, object]]:
    counts: Counter[str] = Counter()
    for seq in alignment.sequences:
        counts.update(ch.upper() for ch in seq if ch.upper() in AA_CANONICAL)
    total = sum(counts.values())
    return [
        {
            "amino_acid": aa,
            "count": counts.get(aa, 0),
            "frequency": counts.get(aa, 0) / total if total else 0.0,
        }
        for aa in sorted(AA_CANONICAL)
    ]


def residue_class_summary(alignment: Alignment) -> list[dict[str, object]]:
    residues = [
        ch.upper()
        for seq in alignment.sequences
        for ch in seq
        if ch.upper() not in MISSING_CHARS and ch.upper() in AA_CANONICAL
    ]
    total = len(residues)
    rows: list[dict[str, object]] = []
    for class_name, members in sorted(RESIDUE_CLASSES.items()):
        count = sum(ch in members for ch in residues)
        rows.append(
            {
                "residue_class": class_name,
                "count": count,
                "fraction": count / total if total else 0.0,
            }
        )
    return rows


def hydrophobicity_summary(alignment: Alignment) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in alignment.records:
        values = [
            HYDROPATHY_KD[ch.upper()] for ch in record.sequence if ch.upper() in HYDROPATHY_KD
        ]
        rows.append(
            {
                "sequence_id": record.id,
                "mean_hydrophobicity_kd": mean(values) if values else None,
                "residue_count": len(values),
            }
        )
    return rows
