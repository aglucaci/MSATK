"""Molecule type detection."""

from __future__ import annotations

from msatk.constants import DNA_ALPHABET, GAP_CHARS, PROTEIN_ALPHABET, RNA_ALPHABET
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
    aliases = {"nucleotide": "dna", "aa": "protein", "amino-acid": "protein", "cds": "codon"}
    return aliases.get(normalized, normalized)
