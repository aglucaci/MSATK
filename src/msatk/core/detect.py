"""Molecule type detection."""

from __future__ import annotations

import numpy as np

from msatk.constants import DNA_ALPHABET, GAP_CHARS, PROTEIN_ALPHABET, RNA_ALPHABET
from msatk.core.genetic_code import genetic_code
from msatk.models import Alignment


def detect_molecule_type(alignment: Alignment) -> str:
    """Return dna, rna, codon, protein, mixed, or unknown."""

    chars = {
        ch.upper()
        for seq in alignment.sequences
        for ch in seq
        if ch.upper() not in GAP_CHARS and ch not in {"?", "*"}
    }
    if not chars:
        return "unknown"
    if chars <= RNA_ALPHABET and "U" in chars and "T" not in chars:
        return "rna"
    if chars <= DNA_ALPHABET:
        return "codon" if alignment.length % 3 == 0 and alignment.length >= 3 else "dna"
    if chars <= PROTEIN_ALPHABET | {"*"}:
        return "protein"
    dna_like = len(chars & DNA_ALPHABET) / max(len(chars), 1)
    protein_like = len(chars & PROTEIN_ALPHABET) / max(len(chars), 1)
    if dna_like > 0.8:
        return "dna"
    if protein_like > 0.8:
        return "protein"
    return "mixed"


def normalize_requested_type(sequence_type: str, alignment: Alignment) -> str:
    if sequence_type == "auto":
        return detect_molecule_type(alignment)
    normalized = sequence_type.lower()
    aliases = {
        "nucleotide": "dna",
        "aa": "protein",
        "amino-acid": "protein",
        "cds": "codon",
        "translated-cds": "translated_cds",
        "translated_cds": "translated_cds",
    }
    return aliases.get(normalized, normalized)


def alignment_detection_summary(
    alignment: Alignment, molecule_type: str, code: int = 1
) -> dict[str, object]:
    """Infer user-facing alignment properties for reports and JSON outputs."""

    length_divisible_by_three = alignment.length % 3 == 0
    nucleotide_like = molecule_type in {"dna", "rna", "codon"}
    stop_rows = _stop_codon_hits(alignment, code) if nucleotide_like else []
    frame_rows = _frame_warnings(alignment) if nucleotide_like else []
    appears_codon_aware = nucleotide_like and length_divisible_by_three and not frame_rows
    return {
        "detected_file_format": alignment.fmt,
        "detected_molecule_type": molecule_type,
        "sequence_lengths_consistent": alignment.is_rectangular,
        "length_divisible_by_three": length_divisible_by_three,
        "appears_codon_aware": appears_codon_aware,
        "stop_codons_exist": bool(stop_rows),
        "stop_codon_count": len(stop_rows),
        "internal_stop_codons_exist": any(bool(row["internal_stop"]) for row in stop_rows),
        "frameshift_warnings_exist": bool(frame_rows),
        "frameshift_warning_count": len(frame_rows),
        "translated_cds_mode": molecule_type == "translated_cds",
        "genetic_code": code if nucleotide_like else None,
    }


def _codon_matrix(alignment: Alignment) -> tuple[np.ndarray, np.ndarray]:
    """Return (uint8 char matrix truncated to whole codons, original lengths)."""

    seqs = alignment.padded_sequences()
    raw = "".join(seqs).upper().replace("U", "T").encode("ascii", errors="replace")
    chars = np.frombuffer(raw, dtype=np.uint8).reshape(len(seqs), alignment.length)
    lengths = np.array([len(r.sequence) for r in alignment.records])
    return chars[:, : (alignment.length // 3) * 3], lengths


def _stop_codon_hits(alignment: Alignment, code: int = 1) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    if alignment.length < 3:
        return hits
    chars, lengths = _codon_matrix(alignment)
    codons = chars.reshape(chars.shape[0], -1, 3)
    is_stop = np.zeros(codons.shape[:2], dtype=bool)
    for stop in sorted(genetic_code(code).stops):
        pattern = np.frombuffer(stop.encode(), dtype=np.uint8)
        is_stop |= (codons == pattern).all(axis=2)
    for row, col in zip(*np.nonzero(is_stop)):
        n_codons = int(lengths[row]) // 3
        if col >= n_codons:
            continue
        hits.append(
            {
                "sequence_id": alignment.records[row].id,
                "codon_index": int(col) + 1,
                "codon": bytes(codons[row, col]).decode(),
                "internal_stop": int(col) + 1 != n_codons,
            }
        )
    return hits


def _frame_warnings(alignment: Alignment) -> list[dict[str, object]]:
    warnings: list[dict[str, object]] = []
    if alignment.length % 3 != 0:
        warnings.append({"scope": "alignment", "reason": "alignment_length_not_divisible_by_three"})
    chars, lengths = _codon_matrix(alignment)
    gap_codes = np.frombuffer("".join(sorted(GAP_CHARS)).encode(), dtype=np.uint8)
    seqs = alignment.padded_sequences()
    full_gap = np.isin(
        np.frombuffer("".join(seqs).encode("ascii", errors="replace"), dtype=np.uint8),
        gap_codes,
    ).reshape(len(seqs), alignment.length)
    ungapped = [int(n) - int(full_gap[i, : int(n)].sum()) for i, n in enumerate(lengths)]
    gaps = full_gap[:, : chars.shape[1]].reshape(chars.shape[0], -1, 3).sum(axis=2)
    partial = (gaps > 0) & (gaps < 3)
    for i, record in enumerate(alignment.records):
        if ungapped[i] % 3 != 0:
            warnings.append(
                {
                    "scope": "sequence",
                    "sequence_id": record.id,
                    "reason": "ungapped_length_not_divisible_by_three",
                    "ungapped_length": ungapped[i],
                }
            )
        n_codons = int(lengths[i]) // 3
        for col in np.nonzero(partial[i, :n_codons])[0]:
            warnings.append(
                {
                    "scope": "sequence",
                    "sequence_id": record.id,
                    "reason": "partial_gap_inside_codon",
                    "codon_index": int(col) + 1,
                }
            )
    return warnings
