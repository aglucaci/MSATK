"""Consensus RNA secondary-structure analysis (Stockholm ``#=GC SS_cons``).

Base pairs are read from WUSS/dot-bracket notation: ``()``, ``<>``, ``[]``, ``{}``
and pseudoknot letter pairs (``A``/``a`` ...). For each consensus pair and each
sequence, the two residues are classified as Watson–Crick (AU, UA, GC, CG), wobble
(GU, UG), non-canonical, or incomplete (gap/ambiguity at either position).
Relative to the most common canonical pair in that column pair, a canonical pair
that differs at both positions is *compensatory* (double substitution, the
classical covariation signal) and one that differs at one position is *consistent*
(e.g. GC→GU).
"""

from __future__ import annotations

from collections import Counter
from typing import Any

from msatk.models import Alignment

_OPEN = {"(": ")", "<": ">", "[": "]", "{": "}"}
_CLOSE = {v: k for k, v in _OPEN.items()}
WATSON_CRICK = {"AU", "UA", "GC", "CG"}
WOBBLE = {"GU", "UG"}


def parse_structure(ss: str) -> list[tuple[int, int]]:
    """Return 1-based ``(i, j)`` base pairs, ``i < j``, from WUSS or dot-bracket."""

    stacks: dict[str, list[int]] = {}
    pairs: list[tuple[int, int]] = []
    for pos, ch in enumerate(ss, start=1):
        if ch in _OPEN:
            stacks.setdefault(ch, []).append(pos)
        elif ch in _CLOSE:
            stack = stacks.get(_CLOSE[ch])
            if not stack:
                raise ValueError(f"Unbalanced '{ch}' at structure position {pos}.")
            pairs.append((stack.pop(), pos))
        elif ch.isalpha() and ch.isupper():
            stacks.setdefault(ch, []).append(pos)
        elif ch.isalpha() and ch.islower():
            stack = stacks.get(ch.upper())
            if not stack:
                raise ValueError(f"Unbalanced pseudoknot '{ch}' at structure position {pos}.")
            pairs.append((stack.pop(), pos))
    left = [k for k, v in stacks.items() if v]
    if left:
        raise ValueError(f"Unclosed structure brackets: {', '.join(left)}.")
    return sorted(pairs)


def _classify(x: str, y: str) -> str:
    pair = (x + y).replace("T", "U")
    if pair in WATSON_CRICK:
        return "watson_crick"
    if pair in WOBBLE:
        return "wobble"
    if x in "ACGUT" and y in "ACGUT":
        return "noncanonical"
    return "incomplete"


def base_pair_rows(alignment: Alignment, ss_cons: str) -> list[dict[str, Any]]:
    pairs = parse_structure(ss_cons)
    seqs = [s.upper() for s in alignment.padded_sequences()]
    rows: list[dict[str, Any]] = []
    for i, j in pairs:
        observed = [(s[i - 1], s[j - 1]) for s in seqs if j <= len(s)]
        classes = Counter(_classify(x, y) for x, y in observed)
        canonical = [
            (x + y).replace("T", "U")
            for x, y in observed
            if _classify(x, y) in {"watson_crick", "wobble"}
        ]
        reference = Counter(canonical).most_common(1)[0][0] if canonical else None
        compensatory = consistent = 0
        if reference:
            for p in canonical:
                diff = (p[0] != reference[0]) + (p[1] != reference[1])
                compensatory += diff == 2
                consistent += diff == 1
        complete = len(observed) - classes["incomplete"]
        rows.append(
            {
                "pair_i": i,
                "pair_j": j,
                "sequences_complete": complete,
                "watson_crick": classes["watson_crick"],
                "wobble": classes["wobble"],
                "noncanonical": classes["noncanonical"],
                "incomplete": classes["incomplete"],
                "canonical_fraction": (classes["watson_crick"] + classes["wobble"]) / complete
                if complete
                else None,
                "dominant_pair": reference,
                "compensatory_sequences": compensatory,
                "consistent_sequences": consistent,
                "covarying": compensatory > 0,
            }
        )
    return rows


def rna_structure_summary(alignment: Alignment, ss_cons: str) -> dict[str, Any]:
    rows = base_pair_rows(alignment, ss_cons)
    fracs = [r["canonical_fraction"] for r in rows if r["canonical_fraction"] is not None]
    paired_sites = {p for r in rows for p in (r["pair_i"], r["pair_j"])}
    return {
        "consensus_base_pairs": len(rows),
        "consensus_paired_sites": len(paired_sites),
        "mean_canonical_pair_fraction": sum(fracs) / len(fracs) if fracs else None,
        "base_pairs_fully_canonical": sum(1 for f in fracs if f == 1.0),
        "base_pairs_with_compensatory_change": sum(1 for r in rows if r["covarying"]),
        "compensatory_pair_observations": sum(r["compensatory_sequences"] for r in rows),
        "consistent_pair_observations": sum(r["consistent_sequences"] for r in rows),
    }
