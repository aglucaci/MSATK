"""Alignment validation helpers."""

from __future__ import annotations

from collections import Counter

from msatk.constants import DNA_ALPHABET, GAP_CHARS, PROTEIN_ALPHABET, RNA_ALPHABET
from msatk.exceptions import AlignmentValidationError
from msatk.models import Alignment

VALID_ALIGNMENT_CHARS = DNA_ALPHABET | RNA_ALPHABET | PROTEIN_ALPHABET | GAP_CHARS | {"?", "*"}


def validate_alignment(alignment: Alignment, mode: str = "permissive") -> list[str]:
    """Validate an alignment and return warnings.

    In ``strict`` mode validation failures raise ``AlignmentValidationError``.
    In ``permissive`` mode the same issues are returned as warnings.
    """

    issues: list[str] = []
    ids = alignment.ids
    duplicates = sorted(seq_id for seq_id, count in Counter(ids).items() if count > 1)
    if duplicates:
        issues.append(f"Duplicate sequence IDs detected: {', '.join(duplicates)}")
    if not alignment.is_rectangular:
        issues.append(
            "Sequence lengths are unequal; MSATK will pad shorter sequences in permissive mode."
        )
    invalid = sorted(
        {ch for seq in alignment.sequences for ch in seq.upper() if ch not in VALID_ALIGNMENT_CHARS}
    )
    if invalid:
        issues.append(f"Invalid alignment character(s) detected: {' '.join(invalid)}")
    if mode == "strict" and issues:
        raise AlignmentValidationError("; ".join(issues))
    return [f"WARNING: {issue}" for issue in issues]
