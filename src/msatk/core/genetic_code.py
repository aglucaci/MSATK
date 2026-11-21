"""NCBI genetic codes (translation tables).

The tables are bundled as ``msatk/data/genetic_codes.json`` (generated from the
NCBI tables distributed with Biopython). For tables 27, 28 and 31, codons that can
be read either as stop or as sense are treated as stops, the conservative choice
for stop-codon and frame QC.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import cache, lru_cache
from importlib import resources
from typing import Any


@dataclass(frozen=True)
class GeneticCode:
    id: int
    name: str
    codons: dict[str, str]  # DNA codon -> one-letter amino acid or '*'
    starts: frozenset[str]
    stops: frozenset[str]

    def translate(self, codon: str) -> str:
        """Translate one codon (U read as T); returns 'X' for non-ACGT codons."""

        return self.codons.get(codon.upper().replace("U", "T"), "X")

    def synonymous_families(self) -> dict[str, list[str]]:
        families: dict[str, list[str]] = {}
        for codon, aa in sorted(self.codons.items()):
            if aa != "*":
                families.setdefault(aa, []).append(codon)
        return families


@lru_cache(maxsize=1)
def _tables() -> dict[str, dict[str, object]]:
    text = resources.files("msatk.data").joinpath("genetic_codes.json").read_text()
    return json.loads(text)["tables"]  # type: ignore[no-any-return]


def available_codes() -> dict[int, str]:
    """Return ``{table id: name}`` for every bundled NCBI genetic code."""

    return {int(k): str(v["name"]) for k, v in _tables().items()}


@cache
def genetic_code(table: int | str = 1) -> GeneticCode:
    """Return an NCBI genetic code by numeric id (default 1, the standard code)."""

    key = str(int(table))
    tables = _tables()
    if key not in tables:
        raise ValueError(
            f"Unknown genetic code {table!r}; available: {', '.join(sorted(tables, key=int))}."
        )
    t: dict[str, Any] = tables[key]
    return GeneticCode(
        id=int(key),
        name=str(t["name"]),
        codons=dict(t["codons"]),
        starts=frozenset(t["starts"]),
        stops=frozenset(t["stops"]),
    )
