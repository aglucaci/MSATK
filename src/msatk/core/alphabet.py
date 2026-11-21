"""Molecule-aware character semantics shared by every MSATK metric.

Every aligned cell is assigned exactly one *state*:

``residue``
    A canonical, unambiguous symbol (``ACGT`` for DNA/codon, ``ACGU`` for RNA,
    the 20 standard amino acids for protein).
``ambiguous``
    A partially informative symbol, e.g. IUPAC ``R``/``Y`` for nucleotides or
    ``B``/``Z``/``J`` for proteins.
``missing``
    Completely uninformative data: ``?`` everywhere, ``N`` for nucleotides and
    ``X`` for proteins. Padding added to unequal-length input is also missing.
``gap``
    An alignment gap (``-`` or ``.``).
``stop``
    A translation stop (``*``) in protein alignments.
``invalid``
    Anything else.

Crucially, the meaning of a character depends on the molecule type: ``N`` is
asparagine in a protein alignment but "any base" in a nucleotide alignment.
"""

from __future__ import annotations

from dataclasses import dataclass, field

GAP_SYMBOLS = frozenset("-.")
PAD_SYMBOL = "?"

NUCLEOTIDE_IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

PROTEIN_AMBIGUITY = {"B": "DN", "Z": "EQ", "J": "IL"}
AA_CANONICAL_ORDER = "ACDEFGHIKLMNPQRSTVWY"

NUCLEOTIDE_TYPES = frozenset({"dna", "rna", "codon"})
PROTEIN_TYPES = frozenset({"protein", "translated_cds"})


@dataclass(frozen=True)
class AlphabetPolicy:
    """How characters are interpreted for one molecule type."""

    molecule_type: str
    canonical: str
    ambiguous: frozenset[str]
    missing: frozenset[str]
    gaps: frozenset[str] = GAP_SYMBOLS
    stop: frozenset[str] = field(default_factory=frozenset)
    synonyms: dict[str, str] = field(default_factory=dict)

    @property
    def canonical_set(self) -> frozenset[str]:
        return frozenset(self.canonical)

    @property
    def alphabet_size(self) -> int:
        return len(self.canonical)

    def normalize(self, ch: str) -> str:
        upper = ch.upper()
        return self.synonyms.get(upper, upper)

    def state(self, ch: str) -> str:
        if ch in self.gaps:
            return "gap"
        symbol = self.normalize(ch)
        if symbol in self.canonical_set:
            return "residue"
        if symbol in self.missing:
            return "missing"
        if symbol in self.ambiguous:
            return "ambiguous"
        if symbol in self.stop:
            return "stop"
        return "invalid"

    def is_residue(self, ch: str) -> bool:
        return self.normalize(ch) in self.canonical_set and ch not in self.gaps

    def describe(self) -> dict[str, object]:
        return {
            "molecule_type": self.molecule_type,
            "canonical": list(self.canonical),
            "ambiguous": sorted(self.ambiguous),
            "missing": sorted(self.missing),
            "gaps": sorted(self.gaps),
            "stop": sorted(self.stop),
            "synonyms": dict(self.synonyms),
            "case_sensitive": False,
        }


DNA_POLICY = AlphabetPolicy(
    molecule_type="dna",
    canonical="ACGT",
    ambiguous=frozenset("RYSWKMBDHV"),
    missing=frozenset({"N", "?"}),
    synonyms={"U": "T"},
)
RNA_POLICY = AlphabetPolicy(
    molecule_type="rna",
    canonical="ACGU",
    ambiguous=frozenset("RYSWKMBDHV"),
    missing=frozenset({"N", "?"}),
    synonyms={"T": "U"},
)
PROTEIN_POLICY = AlphabetPolicy(
    molecule_type="protein",
    canonical=AA_CANONICAL_ORDER,
    ambiguous=frozenset({"B", "Z", "J", "U", "O"}),
    missing=frozenset({"X", "?"}),
    stop=frozenset({"*"}),
)


def policy_for(molecule_type: str) -> AlphabetPolicy:
    """Return the character policy for a molecule type (``mixed``/``unknown`` → DNA-like union)."""

    kind = molecule_type.lower()
    if kind == "rna":
        return RNA_POLICY
    if kind in NUCLEOTIDE_TYPES:
        return AlphabetPolicy(**{**DNA_POLICY.__dict__, "molecule_type": kind})
    if kind in PROTEIN_TYPES:
        return AlphabetPolicy(**{**PROTEIN_POLICY.__dict__, "molecule_type": kind})
    # Mixed/unknown alignments: accept the protein alphabet but treat ?, N and X as missing.
    return AlphabetPolicy(
        molecule_type=kind,
        canonical=AA_CANONICAL_ORDER,
        ambiguous=frozenset({"B", "Z", "J", "U", "O"}),
        missing=frozenset({"?", "X"}),
        stop=frozenset({"*"}),
    )


def is_nucleotide(molecule_type: str) -> bool:
    return molecule_type.lower() in NUCLEOTIDE_TYPES
