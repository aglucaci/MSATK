"""Agreement between alignments of the same sequences.

High conservation does not mean a column is correctly aligned. Alignment
confidence can only be estimated by comparing alignments: against a trusted
reference (benchmarks such as BAliBASE) or against alternative alignments of the
same sequences (different aligners, guide trees or parameters, as in GUIDANCE
or T-Coffee's consistency scores). MSATK keeps these separate from its
descriptive QC scores.

Definitions (Thompson et al. 1999; Edgar 2004):

* An alignment *implies* a residue pair ``(a_p, b_q)`` when residue ``p`` of
  sequence ``a`` and residue ``q`` of sequence ``b`` share a column.
* **SP score** (sum-of-pairs, recall): fraction of the reference's residue pairs
  that the test alignment also implies. **Precision**: fraction of the test's
  pairs implied by the reference. Their harmonic mean is reported as F1.
* **TC score** (total column): fraction of reference columns (with at least two
  residues) reproduced exactly, i.e. the same residues in one test column.
* **Column agreement**: for each column of the test alignment, the fraction of
  its residue pairs implied by the other alignment. Averaged over several
  alternatives it is a per-column confidence score (cf. GUIDANCE).

Sequences are matched by ID, and their residues (gaps removed, case ignored)
must be identical in both alignments.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from msatk.models import Alignment

GAP_CHARS = frozenset("-.")
DEFAULT_MAX_PAIRS = 50000


def _residue_numbers(alignment: Alignment, ids: list[str]) -> tuple[np.ndarray, list[str]]:
    """``(N, L)`` residue number in each cell (``-1`` for gaps) and the ungapped sequences."""

    by_id = {r.id: r.sequence for r in alignment.records}
    seqs = [by_id[i] for i in ids]
    width = max(len(s) for s in seqs)
    out = np.full((len(ids), width), -1, dtype=np.int64)
    ungapped = []
    for r, s in enumerate(seqs):
        is_res = np.array([c not in GAP_CHARS for c in s], dtype=bool)
        out[r, : len(s)][is_res] = np.arange(int(is_res.sum()))
        ungapped.append("".join(c for c in s if c not in GAP_CHARS).upper())
    return out, ungapped


def _check(test: Alignment, other: Alignment) -> list[str]:
    ids_t, ids_o = test.ids, other.ids
    if len(set(ids_t)) != len(ids_t) or len(set(ids_o)) != len(ids_o):
        raise ValueError("Sequence IDs must be unique in both alignments.")
    if set(ids_t) != set(ids_o):
        missing = sorted(set(ids_t) ^ set(ids_o))[:5]
        raise ValueError(f"The alignments contain different sequences, e.g. {missing}.")
    _, a = _residue_numbers(test, ids_t)
    _, b = _residue_numbers(other, ids_t)
    bad = [i for i, x, y in zip(ids_t, a, b) if x != y]
    if bad:
        raise ValueError(
            f"{len(bad)} sequence(s) have different residues in the two alignments, e.g. {bad[:3]}."
        )
    return ids_t


def _partner_maps(res: np.ndarray, i: int, j: int, len_i: int) -> np.ndarray:
    """For each residue of sequence i, the residue of j in the same column (or -1)."""

    both = (res[i] >= 0) & (res[j] >= 0)
    out = np.full(len_i, -1, dtype=np.int64)
    out[res[i][both]] = res[j][both]
    return out


def compare_alignments(
    test: Alignment, reference: Alignment, max_pairs: int = DEFAULT_MAX_PAIRS, seed: int = 0
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    """Score ``test`` against ``reference``: SP, precision, TC and per-column agreement."""

    ids = _check(test, reference)
    rt, ungapped = _residue_numbers(test, ids)
    rr, _ = _residue_numbers(reference, ids)
    n = len(ids)
    lengths = [len(s) for s in ungapped]
    total = n * (n - 1) // 2
    if total <= max_pairs:
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        method = "exact"
    else:
        rng = np.random.default_rng(seed)
        chosen: set[tuple[int, int]] = set()
        while len(chosen) < max_pairs:
            i, j = sorted(rng.choice(n, 2, replace=False).tolist())
            chosen.add((i, j))
        pairs = sorted(chosen)
        method = "sampled"
    col_pairs = np.zeros(rt.shape[1])
    col_agree = np.zeros(rt.shape[1])
    ref_pairs = test_pairs = shared = 0
    seq_shared = np.zeros(n)
    seq_ref = np.zeros(n)
    for i, j in pairs:
        mt = _partner_maps(rt, i, j, lengths[i])
        mr = _partner_maps(rr, i, j, lengths[i])
        in_t, in_r = mt >= 0, mr >= 0
        agree = in_t & in_r & (mt == mr)
        ref_pairs += int(in_r.sum())
        test_pairs += int(in_t.sum())
        shared += int(agree.sum())
        seq_shared[[i, j]] += agree.sum()
        seq_ref[[i, j]] += in_r.sum()
        # Column agreement for the test alignment's columns.
        both = (rt[i] >= 0) & (rt[j] >= 0)
        cols = np.flatnonzero(both)
        col_pairs[cols] += 1
        col_agree[cols] += agree[rt[i][cols]]
    sp = shared / ref_pairs if ref_pairs else None
    precision = shared / test_pairs if test_pairs else None
    f1 = 2 * sp * precision / (sp + precision) if sp and precision else None

    # TC: reference columns reproduced exactly (all sequences, not sampled).
    def signatures(res: np.ndarray) -> list[frozenset[tuple[int, int]]]:
        out = []
        for c in range(res.shape[1]):
            rows = np.flatnonzero(res[:, c] >= 0)
            if rows.size >= 2:
                out.append(frozenset(zip(rows.tolist(), res[rows, c].tolist())))
        return out

    ref_cols = signatures(rr)
    test_cols = set(signatures(rt))
    tc = sum(c in test_cols for c in ref_cols) / len(ref_cols) if ref_cols else None
    column_rows = [
        {
            "site_index": c + 1,
            "residue_pairs": int(col_pairs[c]),
            "agreement": float(col_agree[c] / col_pairs[c]) if col_pairs[c] else None,
        }
        for c in range(rt.shape[1])
    ]
    seq_rows = [
        {
            "sequence_id": seq_id,
            "pairs_recovered": float(seq_shared[k] / seq_ref[k]) if seq_ref[k] else None,
        }
        for k, seq_id in enumerate(ids)
    ]
    agreements = np.array([r["agreement"] for r in column_rows if r["agreement"] is not None])
    summary = {
        "pair_method": method,
        "sequence_pairs": len(pairs),
        "reference_residue_pairs": ref_pairs,
        "test_residue_pairs": test_pairs,
        "shared_residue_pairs": shared,
        "sp_score": sp,
        "precision": precision,
        "f1": f1,
        "tc_score": tc,
        "reference_columns": len(ref_cols),
        "mean_column_agreement": float(agreements.mean()) if agreements.size else None,
        "columns_below_0_5_agreement": int((agreements < 0.5).sum()),
    }
    return summary, {"column_agreement": column_rows, "sequence_agreement": seq_rows}


def ensemble_agreement(
    test: Alignment,
    alternatives: list[Alignment],
    max_pairs: int = DEFAULT_MAX_PAIRS,
    seed: int = 0,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    """Mean column agreement of ``test`` with each alternative (a per-column confidence)."""

    if not alternatives:
        raise ValueError("Give at least one alternative alignment.")
    results = [compare_alignments(test, alt, max_pairs, seed) for alt in alternatives]
    width = len(results[0][1]["column_agreement"])
    rows = []
    for c in range(width):
        vals = [r[1]["column_agreement"][c]["agreement"] for r in results]
        vals = [v for v in vals if v is not None]
        rows.append(
            {
                "site_index": c + 1,
                "confidence": float(np.mean(vals)) if vals else None,
                "min_agreement": float(np.min(vals)) if vals else None,
            }
        )
    conf = np.array([r["confidence"] for r in rows if r["confidence"] is not None])
    summary = {
        "alternatives": len(alternatives),
        "mean_sp_vs_alternatives": float(np.mean([r[0]["sp_score"] or 0 for r in results])),
        "mean_tc_vs_alternatives": float(np.mean([r[0]["tc_score"] or 0 for r in results])),
        "mean_column_confidence": float(conf.mean()) if conf.size else None,
        "columns_below_0_5_confidence": int((conf < 0.5).sum()),
        "columns_fully_consistent": int((conf == 1.0).sum()),
    }
    return summary, {"column_confidence": rows}
