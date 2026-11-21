"""Analyses that apply only to some sequence alignments.

Each analysis runs only when it is biologically and mathematically valid for the
input: pairwise statistics for exactly two sequences, substitution-matrix scores
for proteins, base-pair analysis when a consensus structure is present, profile
state statistics when match/insert states are known, and per-partition statistics
when partitions are declared. ``applicability`` reports why each was or was not
run, so nothing is silently skipped.
"""

from __future__ import annotations

from typing import Any

from msatk.core.alphabet import is_nucleotide
from msatk.io.readers import Partition, SequenceAnnotations
from msatk.stats.alignment import StatisticsBundle


def specialized_analyses(
    bundle: StatisticsBundle,
    molecule_type: str,
    annotations: SequenceAnnotations | None = None,
    partitions: list[Partition] | None = None,
) -> tuple[dict[str, dict[str, Any]], dict[str, list[dict[str, Any]]], dict[str, str]]:
    """Return ``(summaries, tables, applicability)`` keyed by analysis name."""

    ann = annotations or SequenceAnnotations()
    alignment = bundle.encoded.alignment
    summaries: dict[str, dict[str, Any]] = {}
    tables: dict[str, list[dict[str, Any]]] = {}
    why: dict[str, str] = {}

    if alignment.n_sequences == 2:
        from msatk.stats.pairwise_alignment import pairwise_alignment_stats

        a, b = alignment.padded_sequences()
        summaries["pairwise_alignment"] = pairwise_alignment_stats(
            a, b, molecule_type, method=ann.method, reported=ann.pairwise or None
        )
        why["pairwise_alignment"] = "run: exactly two sequences"
    else:
        why["pairwise_alignment"] = f"skipped: {alignment.n_sequences} sequences (needs exactly 2)"

    if molecule_type in {"protein", "translated_cds"}:
        from msatk.stats.substitution import column_substitution_scores, substitution_summary

        summaries["substitution"] = substitution_summary(bundle.encoded)
        tables["site_substitution_scores"] = column_substitution_scores(bundle.encoded)
        why["substitution"] = "run: protein alignment (BLOSUM62)"
    else:
        why["substitution"] = "skipped: substitution matrices apply to proteins"

    if ann.ss_cons:
        if is_nucleotide(molecule_type):
            from msatk.stats.rna_structure import base_pair_rows, rna_structure_summary

            try:
                summaries["rna_structure"] = rna_structure_summary(alignment, ann.ss_cons)
                tables["rna_base_pairs"] = base_pair_rows(alignment, ann.ss_cons)
                why["rna_structure"] = "run: #=GC SS_cons present"
            except ValueError as exc:
                why["rna_structure"] = f"skipped: invalid SS_cons ({exc})"
        else:
            why["rna_structure"] = "skipped: SS_cons present but alignment is not nucleotide"
    else:
        why["rna_structure"] = "skipped: no consensus structure (Stockholm #=GC SS_cons)"

    if ann.insertions is not None:
        from msatk.stats.profile_states import (
            profile_rows_from_insertions,
            profile_summary_from_insertions,
        )

        summaries["profile_states"] = profile_summary_from_insertions(alignment, ann.insertions)
        tables["profile_states"] = profile_rows_from_insertions(alignment, ann.insertions)
        why["profile_states"] = "run: A2M/A3M insert states"
    elif ann.rf:
        from msatk.stats.profile_states import profile_summary_from_rf

        summaries["profile_states"] = profile_summary_from_rf(alignment, ann.rf)
        why["profile_states"] = "run: Stockholm #=GC RF match columns"
    else:
        why["profile_states"] = "skipped: no match/insert state annotation (A2M/A3M or #=GC RF)"

    parts = list(partitions or []) or ann.partitions
    if parts:
        from msatk.stats.partitions import partition_rows, partition_summary

        rows = partition_rows(alignment, molecule_type, parts)
        summaries["partitions"] = partition_summary(parts, rows, alignment.length)
        tables["partitions"] = rows
        why["partitions"] = f"run: {len(parts)} partitions declared"
    else:
        why["partitions"] = "skipped: no partitions (NEXUS charsets or --partitions file)"
    return summaries, tables, why
