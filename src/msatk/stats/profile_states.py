"""Profile (match/insert state) statistics for A2M/A3M and Stockholm ``RF`` input.

A2M/A3M: uppercase residues and ``-`` are match states; lowercase residues are
insertions and are removed from the analysed matrix, so insert lengths come from
the parser. Stockholm: the ``#=GC RF`` reference line marks match columns
(``x``/residue letters) and insert columns (``.``/``-``).
"""

from __future__ import annotations

from typing import Any

from msatk.models import Alignment


def match_columns_from_rf(rf: str) -> list[bool]:
    return [ch not in ".-~" for ch in rf]


def profile_rows_from_insertions(
    alignment: Alignment, insertions: list[list[int]]
) -> list[dict[str, Any]]:
    """Per-sequence match-state coverage and insertion totals (A2M/A3M)."""

    rows = []
    for rec, ins in zip(alignment.records, insertions):
        match = len(rec.sequence)
        occupied = sum(1 for ch in rec.sequence if ch not in "-.")
        rows.append(
            {
                "sequence_id": rec.id,
                "match_states": match,
                "match_states_occupied": occupied,
                "match_state_coverage": occupied / match if match else None,
                "deletions": match - occupied,
                "inserted_residues": sum(ins),
                "insertion_events": sum(1 for v in ins if v),
            }
        )
    return rows


def profile_summary_from_insertions(
    alignment: Alignment, insertions: list[list[int]]
) -> dict[str, Any]:
    rows = profile_rows_from_insertions(alignment, insertions)
    n = len(rows)
    width = alignment.length
    per_column = [sum(ins[k] > 0 for ins in insertions) for k in range(width + 1)] if n else []
    return {
        "profile_source": "a2m/a3m insert states",
        "match_states": width,
        "mean_match_state_coverage": sum(r["match_state_coverage"] or 0 for r in rows) / n
        if n
        else None,
        "total_inserted_residues": sum(r["inserted_residues"] for r in rows),
        "sequences_with_insertions": sum(1 for r in rows if r["inserted_residues"]),
        "insert_positions_used": sum(1 for v in per_column if v),
        "max_sequences_inserting_at_one_position": max(per_column) if per_column else 0,
    }


def profile_summary_from_rf(alignment: Alignment, rf: str) -> dict[str, Any]:
    match = match_columns_from_rf(rf)
    seqs = alignment.padded_sequences()
    n_match = sum(match)
    cov = []
    inserted = []
    for s in seqs:
        occ = sum(1 for k, ch in enumerate(s) if k < len(match) and match[k] and ch not in "-.")
        cov.append(occ / n_match if n_match else 0.0)
        inserted.append(
            sum(1 for k, ch in enumerate(s) if k < len(match) and not match[k] and ch not in "-.")
        )
    return {
        "profile_source": "Stockholm #=GC RF",
        "match_states": n_match,
        "insert_columns": len(match) - n_match,
        "mean_match_state_coverage": sum(cov) / len(cov) if cov else None,
        "total_inserted_residues": sum(inserted),
        "sequences_with_insertions": sum(1 for v in inserted if v),
    }
