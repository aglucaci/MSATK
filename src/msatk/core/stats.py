"""Alignment statistics."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable
from statistics import mean, median

from msatk.constants import DNA_BASES, GAP_CHARS, MISSING_CHARS
from msatk.models import Alignment


def shannon_entropy(values: Iterable[str]) -> float:
    counts = Counter(v.upper() for v in values if v.upper() not in MISSING_CHARS)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def gc_content(sequence: str) -> float | None:
    chars = [ch.upper().replace("U", "T") for ch in sequence if ch.upper() not in GAP_CHARS]
    bases = [ch for ch in chars if ch in DNA_BASES]
    if not bases:
        return None
    return (bases.count("G") + bases.count("C")) / len(bases)


def pairwise_identity(seq_a: str, seq_b: str) -> float:
    matches = 0
    comparable = 0
    for a, b in zip(seq_a, seq_b):
        if a.upper() in MISSING_CHARS or b.upper() in MISSING_CHARS:
            continue
        comparable += 1
        matches += int(a.upper() == b.upper())
    return matches / comparable if comparable else 0.0


def identity_matrix(alignment: Alignment) -> list[list[float]]:
    seqs = alignment.padded_sequences()
    return [[pairwise_identity(a, b) for b in seqs] for a in seqs]


def distance_matrix(alignment: Alignment) -> list[list[float]]:
    return [[1.0 - value for value in row] for row in identity_matrix(alignment)]


def per_sequence_stats(alignment: Alignment) -> list[dict[str, object]]:
    seqs = alignment.padded_sequences()
    ids = alignment.ids
    matrix = identity_matrix(alignment)
    rows: list[dict[str, object]] = []
    for idx, (seq_id, seq) in enumerate(zip(ids, seqs)):
        raw_length = len(seq.rstrip("-"))
        gap_count = sum(ch in GAP_CHARS for ch in seq)
        ambiguous_count = sum(ch.upper() in MISSING_CHARS and ch not in GAP_CHARS for ch in seq)
        ungapped = "".join(ch for ch in seq if ch not in GAP_CHARS)
        peer_identities = [value for j, value in enumerate(matrix[idx]) if j != idx]
        mean_identity = mean(peer_identities) if peer_identities else 1.0
        rows.append(
            {
                "sequence_id": seq_id,
                "raw_length": raw_length,
                "ungapped_length": len(ungapped),
                "gap_count": gap_count,
                "gap_fraction": gap_count / len(seq) if seq else 0.0,
                "ambiguous_count": ambiguous_count,
                "ambiguous_fraction": ambiguous_count / len(seq) if seq else 0.0,
                "gc_content": gc_content(seq),
                "unique_residue_count": len(
                    {ch.upper() for ch in ungapped if ch.upper() not in MISSING_CHARS}
                ),
                "mean_identity_to_others": mean_identity,
                "outlier_score": 1.0 - mean_identity,
            }
        )
    return rows


def per_site_stats(alignment: Alignment) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index, column in enumerate(alignment.columns(), start=1):
        counts = Counter(ch.upper() for ch in column if ch.upper() not in MISSING_CHARS)
        total = len(column)
        observed = sum(counts.values())
        gap_count = sum(ch in GAP_CHARS for ch in column)
        dominant, dominant_count = counts.most_common(1)[0] if counts else ("", 0)
        informative_states = sum(1 for count in counts.values() if count >= 2)
        rows.append(
            {
                "site_index": index,
                "dominant_character": dominant,
                "dominant_character_frequency": dominant_count / observed if observed else 0.0,
                "gap_fraction": gap_count / total if total else 0.0,
                "missing_fraction": sum(ch.upper() in MISSING_CHARS for ch in column) / total
                if total
                else 0.0,
                "entropy": shannon_entropy(column),
                "conservation_score": dominant_count / observed if observed else 0.0,
                "observed_state_count": len(counts),
                "variable": len(counts) > 1,
                "parsimony_informative": informative_states >= 2,
                "singleton": len(counts) > 1 and any(count == 1 for count in counts.values()),
                "codon_position": ((index - 1) % 3) + 1,
            }
        )
    return rows


def alignment_summary(alignment: Alignment, molecule_type: str) -> dict[str, object]:
    seq_rows = per_sequence_stats(alignment)
    site_rows = per_site_stats(alignment)
    identities = [
        value
        for i, row in enumerate(identity_matrix(alignment))
        for j, value in enumerate(row)
        if i < j
    ]
    gaps = sum(row["gap_count"] for row in seq_rows)
    cells = alignment.n_sequences * alignment.length
    gc_values = [row["gc_content"] for row in seq_rows if row["gc_content"] is not None]
    return {
        "tool": "MSATK",
        "number_of_sequences": alignment.n_sequences,
        "alignment_length": alignment.length,
        "rectangular_alignment": alignment.is_rectangular,
        "detected_format": alignment.fmt,
        "molecule_type": molecule_type,
        "mean_sequence_length": mean(row["ungapped_length"] for row in seq_rows),
        "median_sequence_length": median(row["ungapped_length"] for row in seq_rows),
        "total_gaps": gaps,
        "gap_fraction": gaps / cells if cells else 0.0,
        "conserved_sites": sum(not row["variable"] for row in site_rows),
        "variable_sites": sum(bool(row["variable"]) for row in site_rows),
        "parsimony_informative_sites": sum(bool(row["parsimony_informative"]) for row in site_rows),
        "singleton_sites": sum(bool(row["singleton"]) for row in site_rows),
        "mean_pairwise_identity": mean(identities) if identities else 1.0,
        "mean_pairwise_distance": 1.0 - mean(identities) if identities else 0.0,
        "entropy_mean": mean(row["entropy"] for row in site_rows) if site_rows else 0.0,
        "entropy_median": median(row["entropy"] for row in site_rows) if site_rows else 0.0,
        "gc_content": mean(gc_values) if gc_values else None,
    }


def composition_summary(alignment: Alignment) -> list[dict[str, object]]:
    counts: Counter[str] = Counter()
    for seq in alignment.sequences:
        counts.update(ch.upper() for ch in seq if ch.upper() not in MISSING_CHARS)
    total = sum(counts.values())
    return [
        {"character": char, "count": count, "fraction": count / total if total else 0.0}
        for char, count in sorted(counts.items())
    ]
