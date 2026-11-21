"""Codon usage, RSCU, and translation summaries."""

from __future__ import annotations

from collections import Counter, defaultdict
from statistics import mean
from typing import Any

from msatk.constants import GAP_CHARS
from msatk.core.genetic_code import genetic_code
from msatk.models import Alignment


def iter_codons(sequence: str) -> list[tuple[int, str]]:
    clean = sequence.upper().replace("U", "T")
    return [(i // 3 + 1, clean[i : i + 3]) for i in range(0, len(clean) - 2, 3)]


def codon_usage(alignment: Alignment, code: int = 1) -> list[dict[str, Any]]:
    table = genetic_code(code).codons
    counts: Counter[str] = Counter()
    for seq in alignment.padded_sequences():
        for _, codon in iter_codons(seq):
            if len(codon) == 3 and not any(ch in GAP_CHARS or ch == "N" for ch in codon):
                counts[codon] += 1
    total = sum(counts.values())
    return [
        {
            "codon": codon,
            "amino_acid": table.get(codon, "X"),
            "count": counts.get(codon, 0),
            "frequency": counts.get(codon, 0) / total if total else 0.0,
        }
        for codon in sorted(table)
    ]


def rscu(alignment: Alignment, code: int = 1) -> list[dict[str, Any]]:
    usage = {row["codon"]: int(row["count"]) for row in codon_usage(alignment, code)}
    by_aa: dict[str, list[str]] = defaultdict(list)
    for aa, codons in genetic_code(code).synonymous_families().items():
        by_aa[aa].extend(codons)
    rows: list[dict[str, Any]] = []
    for aa, codons in sorted(by_aa.items()):
        family_total = sum(usage.get(codon, 0) for codon in codons)
        expected = family_total / len(codons) if codons else 0
        for codon in sorted(codons):
            count = usage.get(codon, 0)
            rows.append(
                {
                    "amino_acid": aa,
                    "codon": codon,
                    "count": count,
                    "synonymous_family_size": len(codons),
                    "rscu": count / expected if expected else None,
                }
            )
    return rows


def amino_acid_from_codons(alignment: Alignment, code: int = 1) -> list[dict[str, Any]]:
    table = genetic_code(code).codons
    counts: Counter[str] = Counter()
    for seq in alignment.padded_sequences():
        for _, codon in iter_codons(seq):
            if len(codon) == 3 and not any(ch in GAP_CHARS or ch == "N" for ch in codon):
                counts[table.get(codon, "X")] += 1
    total = sum(counts.values())
    return [
        {"amino_acid": aa, "count": count, "frequency": count / total if total else 0.0}
        for aa, count in sorted(counts.items())
    ]


def stop_codon_report(alignment: Alignment, code: int = 1) -> list[dict[str, Any]]:
    stops = genetic_code(code).stops
    rows: list[dict[str, Any]] = []
    for record in alignment.records:
        codons = iter_codons(record.sequence)
        for index, codon in codons:
            if codon.replace("U", "T") in stops:
                rows.append(
                    {
                        "sequence_id": record.id,
                        "codon_index": index,
                        "codon": codon,
                        "internal_stop": index != len(codons),
                    }
                )
    return rows


def codon_position_gc(alignment: Alignment) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for position in (1, 2, 3):
        chars = [
            seq[i].upper().replace("U", "T")
            for seq in alignment.padded_sequences()
            for i in range(position - 1, len(seq), 3)
            if seq[i].upper() in {"A", "C", "G", "T", "U"}
        ]
        rows.append(
            {
                "codon_position": position,
                "gc_content": (chars.count("G") + chars.count("C")) / len(chars) if chars else None,
                "observed_bases": len(chars),
            }
        )
    return rows


def codon_position_entropy(per_site: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for position in (1, 2, 3):
        values = [
            float(row["entropy"])
            for row in per_site
            if row.get("codon_position") == position and row.get("entropy") is not None
        ]
        rows.append(
            {
                "codon_position": position,
                "mean_entropy": mean(values) if values else None,
                "observed_sites": len(values),
            }
        )
    return rows
