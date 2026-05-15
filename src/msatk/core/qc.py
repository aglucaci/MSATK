"""Quality-control warnings and flags."""

from __future__ import annotations

from typing import Any

from msatk.models import Alignment


def flagged_sequences(
    per_sequence: list[dict[str, Any]], max_gap_seq: float = 0.3
) -> list[dict[str, Any]]:
    return [row for row in per_sequence if float(row.get("gap_fraction", 0.0)) > max_gap_seq]


def flagged_sites(
    per_site: list[dict[str, Any]], max_gap_site: float = 0.5
) -> list[dict[str, Any]]:
    return [row for row in per_site if float(row.get("gap_fraction", 0.0)) > max_gap_site]


def qc_warnings(
    alignment: Alignment,
    molecule_type: str,
    per_sequence: list[dict[str, Any]],
    per_site: list[dict[str, Any]],
    max_gap_seq: float = 0.3,
    max_gap_site: float = 0.5,
) -> list[str]:
    warnings: list[str] = []
    seq_flags = flagged_sequences(per_sequence, max_gap_seq)
    site_flags = flagged_sites(per_site, max_gap_site)
    if not alignment.is_rectangular:
        warnings.append(
            "WARNING: Sequence lengths are inconsistent; MSATK padded shorter sequences for site-level summaries."
        )
    if seq_flags:
        warnings.append(f"WARNING: {len(seq_flags)} sequences have >{max_gap_seq:.0%} gaps.")
    if site_flags:
        warnings.append(
            f"WARNING: {len(site_flags)} alignment columns have >{max_gap_site:.0%} missing data."
        )
    if molecule_type == "codon" and alignment.length % 3 != 0:
        warnings.append(
            "WARNING: Alignment length is not divisible by 3; codon-aware analysis may be invalid."
        )
    if molecule_type == "codon":
        for record in alignment.records:
            seq = record.sequence.upper().replace("U", "T")
            codons = [seq[i : i + 3] for i in range(0, len(seq) - 2, 3)]
            internal_stops = [codon for codon in codons[:-1] if codon in {"TAA", "TAG", "TGA"}]
            if internal_stops:
                warnings.append(
                    f"WARNING: Sequence {record.id} contains {len(internal_stops)} internal stop codon(s)."
                )
    outliers = [row for row in per_sequence if float(row.get("outlier_score", 0.0)) > 0.5]
    if outliers:
        warnings.append(
            f"WARNING: {len(outliers)} sequences have unusually low mean identity to the alignment."
        )
    if not warnings:
        warnings.append("OK: No major MSATK QC warnings triggered by current thresholds.")
    return warnings
