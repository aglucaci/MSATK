"""Statistics for a two-sequence (pairwise) alignment.

These are the quantities reported by pairwise aligners: identity, similarity,
mismatches, gap openings and extensions, terminal versus internal gaps, coverage
of each sequence and a recomputed alignment score under an explicit scheme. The
alignment method (global, local, semiglobal) is recorded when the input format
states it (EMBOSS) and is never inferred from the final alignment.
"""

from __future__ import annotations

from typing import Any

from msatk.core.alphabet import is_nucleotide, policy_for
from msatk.stats.substitution import blosum62_score

#: EMBOSS defaults: EDNAFULL-style +5/-4 for nucleotides, BLOSUM62 for proteins,
#: gap open 10 (first gap position) and extend 0.5, end gaps not penalised.
DEFAULT_SCORING = {"match": 5.0, "mismatch": -4.0, "gap_open": 10.0, "gap_extend": 0.5}


def _gap_runs(seq: str, gaps: frozenset[str]) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    start = None
    for i, ch in enumerate(seq):
        if ch in gaps and start is None:
            start = i
        elif ch not in gaps and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(seq)))
    return runs


def pairwise_alignment_stats(
    seq_a: str,
    seq_b: str,
    molecule_type: str,
    method: str = "unknown",
    scoring: dict[str, float] | None = None,
    reported: dict[str, Any] | None = None,
) -> dict[str, Any]:
    policy = policy_for(molecule_type)
    a, b = seq_a.upper(), seq_b.upper()
    width = max(len(a), len(b))
    a, b = a.ljust(width, "?"), b.ljust(width, "?")
    gaps = policy.gaps
    canon = policy.canonical_set
    nucleotide = is_nucleotide(molecule_type)
    sc = {**DEFAULT_SCORING, **(scoring or {})}

    identical = mismatches = similar = uncertain = 0
    for x, y in zip(a, b):
        if x in gaps or y in gaps:
            continue
        if x in canon and y in canon:
            if policy.normalize(x) == policy.normalize(y):
                identical += 1
            else:
                mismatches += 1
                if not nucleotide and (blosum62_score(x, y) or 0) > 0:
                    similar += 1
        else:
            uncertain += 1
    compared = identical + mismatches
    out: dict[str, Any] = {
        "alignment_method": method,
        "alignment_columns": width,
        "aligned_pairs": compared + uncertain,
        "identical_pairs": identical,
        "mismatched_pairs": mismatches,
        "uncertain_pairs": uncertain,
        "identity": identical / compared if compared else None,
        "identity_over_columns": identical / width if width else None,
        "similarity": (identical + similar) / compared if compared and not nucleotide else None,
    }
    gap_open = gap_ext = terminal = internal = 0
    for k, seq in ((1, a), (2, b)):
        runs = _gap_runs(seq, gaps)
        residues = sum(1 for ch in seq if ch not in gaps and ch != "?")
        term = [r for r in runs if r[0] == 0 or r[1] == width]
        intern = [r for r in runs if r not in term]
        gap_open += len(intern)
        gap_ext += sum(e - s - 1 for s, e in intern)
        terminal += sum(e - s for s, e in term)
        internal += sum(e - s for s, e in intern)
        other = b if k == 1 else a
        paired = sum(1 for x, y in zip(seq, other) if x not in gaps and x != "?" and y not in gaps)
        out[f"length_{k}"] = residues
        out[f"coverage_{k}"] = paired / residues if residues else None
    out.update(
        {
            "gap_openings": gap_open,
            "gap_extensions": gap_ext,
            "internal_gap_positions": internal,
            "terminal_gap_positions": terminal,
        }
    )
    # Score with end gaps free (EMBOSS default); uncertain pairs score 0.
    score = 0.0
    for x, y in zip(a, b):
        if x in gaps or y in gaps:
            continue
        if nucleotide:
            if x in canon and y in canon:
                score += (
                    sc["match"] if policy.normalize(x) == policy.normalize(y) else sc["mismatch"]
                )
        else:
            score += blosum62_score(x, y) or 0
    for _, (s, e) in ((None, r) for seq in (a, b) for r in _gap_runs(seq, gaps)):
        if s == 0 or e == width:
            continue
        score -= sc["gap_open"] + sc["gap_extend"] * (e - s - 1)
    out["score_recomputed"] = score
    out["scoring_scheme"] = (
        f"{'match ' + str(sc['match']) + '/mismatch ' + str(sc['mismatch']) if nucleotide else 'BLOSUM62'}"
        f", gap open {sc['gap_open']}, extend {sc['gap_extend']}, end gaps free"
    )
    if reported:
        out["score_reported"] = reported.get("score")
        out["program_reported"] = reported.get("program")
    return out
