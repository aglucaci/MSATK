"""Molecule type detection."""

from __future__ import annotations

from msatk.constants import DNA_ALPHABET, GAP_CHARS, PROTEIN_ALPHABET, RNA_ALPHABET, STOP_CODONS
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


def alignment_detection_summary(alignment: Alignment, molecule_type: str) -> dict[str, object]:
    """Infer user-facing alignment properties for reports and JSON outputs."""

    length_divisible_by_three = alignment.length % 3 == 0
    nucleotide_like = molecule_type in {"dna", "rna", "codon"}
    stop_rows = _stop_codon_hits(alignment) if nucleotide_like else []
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
    }


def _stop_codon_hits(alignment: Alignment) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    if alignment.length < 3:
        return hits
    for record in alignment.records:
        seq = record.sequence.upper().replace("U", "T")
        codons = [seq[i : i + 3] for i in range(0, len(seq) - 2, 3)]
        for index, codon in enumerate(codons, start=1):
            if codon in STOP_CODONS:
                hits.append(
                    {
                        "sequence_id": record.id,
                        "codon_index": index,
                        "codon": codon,
                        "internal_stop": index != len(codons),
                    }
                )
    return hits


def _frame_warnings(alignment: Alignment) -> list[dict[str, object]]:
    warnings: list[dict[str, object]] = []
    if alignment.length % 3 != 0:
        warnings.append({"scope": "alignment", "reason": "alignment_length_not_divisible_by_three"})
    for record in alignment.records:
        ungapped = "".join(ch for ch in record.sequence if ch not in GAP_CHARS)
        if len(ungapped) % 3 != 0:
            warnings.append(
                {
                    "scope": "sequence",
                    "sequence_id": record.id,
                    "reason": "ungapped_length_not_divisible_by_three",
                    "ungapped_length": len(ungapped),
                }
            )
        seq = record.sequence.upper()
        for index in range(0, max(len(seq) - 2, 0), 3):
            codon = seq[index : index + 3]
            if any(ch in GAP_CHARS for ch in codon) and not all(ch in GAP_CHARS for ch in codon):
                warnings.append(
                    {
                        "scope": "sequence",
                        "sequence_id": record.id,
                        "reason": "partial_gap_inside_codon",
                        "codon_index": index // 3 + 1,
                    }
                )
    return warnings
