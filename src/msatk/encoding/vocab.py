"""Fixed, documented vocabularies for categorical and one-hot alignment encodings."""

from __future__ import annotations

from dataclasses import dataclass

from msatk.core.alphabet import (
    AA_CANONICAL_ORDER,
    NUCLEOTIDE_IUPAC,
    PROTEIN_AMBIGUITY,
    is_nucleotide,
    policy_for,
)


@dataclass(frozen=True)
class Vocabulary:
    """An ordered symbol vocabulary. Index ``len(symbols)`` is reserved for 'other'."""

    name: str
    symbols: tuple[str, ...]

    @property
    def other_index(self) -> int:
        return len(self.symbols)

    @property
    def size(self) -> int:
        return len(self.symbols) + 1

    def index(self) -> dict[str, int]:
        return {s: i for i, s in enumerate(self.symbols)}

    def describe(self) -> dict[str, object]:
        return {
            "name": self.name,
            "symbols": list(self.symbols),
            "other_index": self.other_index,
            "case_folded": True,
        }


NUCLEOTIDE_FULL = Vocabulary("nucleotide_iupac", tuple("ACGTURYSWKMBDHVN") + ("-", ".", "?", "*"))
PROTEIN_FULL = Vocabulary(
    "protein_full", tuple(AA_CANONICAL_ORDER) + tuple("BZJUOX") + ("*", "-", ".", "?")
)


def full_vocabulary(molecule_type: str) -> Vocabulary:
    return NUCLEOTIDE_FULL if is_nucleotide(molecule_type) else PROTEIN_FULL


def ambiguity_expansion(molecule_type: str) -> dict[str, str]:
    """Canonical residues compatible with each ambiguity symbol (``N``/``X`` → all residues)."""

    if is_nucleotide(molecule_type):
        policy = policy_for(molecule_type)
        canon = policy.canonical
        expanded = {}
        for code, bases in NUCLEOTIDE_IUPAC.items():
            if code in canon or code in {"U", "T"}:
                continue
            expanded[code] = "".join(policy.normalize(b) for b in bases)
        return expanded
    return {**PROTEIN_AMBIGUITY, "X": AA_CANONICAL_ORDER}
