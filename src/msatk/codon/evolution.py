"""Codon-level evolutionary statistics for coding alignments.

* **Pairwise dN, dS and dN/dS** by the Nei–Gojobori (1986) method with the
  Jukes–Cantor correction, under the selected genetic code. For codons that
  differ at two or three positions, every mutational pathway that avoids stop
  codons is weighted equally, as in the original method. (Biopython's ``NG86``
  also counts pathways through stop codons, scoring a stop-to-stop step as
  synonymous; with those two rules MSATK reproduces it to 1e-10, which the test
  suite checks.) Codons containing gaps, ambiguity codes or stops
  are excluded pairwise. This is a counting estimate between two sequences, not
  a phylogenetic test of selection: use codon models (HyPhy, PAML) for that.
  dN/dS is reported as undefined when dS is 0, saturated (p_S >= 0.75) or
  based on fewer than ``MIN_CODONS`` codons.
* **Effective number of codons, ENC** (Wright 1990), from all complete codons,
  grouping amino acids by synonymous family size under the selected code;
  missing family classes are interpolated as Wright describes (``F3`` from
  ``F2`` and ``F4``), and ENC is capped at 61.
* **Synonymous and non-synonymous polymorphic codon sites**: variable codon
  columns whose observed codons all encode the same amino acid, versus those
  with more than one amino acid.
* **Four-fold degenerate sites**: third positions of codon columns whose
  consensus codon is four-fold degenerate there; reports their nucleotide
  diversity (mean pairwise difference per site).
"""

from __future__ import annotations

import itertools
import math
from functools import cache
from typing import Any

import numpy as np

from msatk.core.genetic_code import GeneticCode, genetic_code

BASES = "TCAG"
CODONS = ["".join(c) for c in itertools.product(BASES, repeat=3)]
CODON_INDEX = {c: i for i, c in enumerate(CODONS)}
MIN_CODONS = 10
DEFAULT_MAX_PAIRS = 20000


@cache
def _tables(code_id: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-codon synonymous sites, and per codon-pair (Sd, Nd, valid) tables."""

    code = genetic_code(code_id)
    aa = [code.translate(c) for c in CODONS]
    stop = np.array([a == "*" for a in aa])
    syn_sites = np.zeros(64)
    for i, c in enumerate(CODONS):
        if stop[i]:
            continue
        s = 0.0
        for pos in range(3):
            alts = [c[:pos] + b + c[pos + 1 :] for b in BASES if b != c[pos]]
            nonstop = [a for a in alts if code.translate(a) != "*"]
            # Nei–Gojobori: fraction of the three possible changes that are synonymous
            s += sum(code.translate(a) == aa[i] for a in nonstop) / 3.0
        syn_sites[i] = s
    sd = np.zeros((64, 64))
    nd = np.zeros((64, 64))
    valid = np.zeros((64, 64), dtype=bool)
    for i, j in itertools.product(range(64), repeat=2):
        if stop[i] or stop[j]:
            continue
        a, b = CODONS[i], CODONS[j]
        diff = [p for p in range(3) if a[p] != b[p]]
        if not diff:
            valid[i, j] = True
            continue
        paths_s, paths_n, paths = 0.0, 0.0, 0
        for order in itertools.permutations(diff):
            cur = a
            ps = pn = 0
            ok = True
            for p in order:
                nxt = cur[:p] + b[p] + cur[p + 1 :]
                if code.translate(nxt) == "*":
                    ok = False
                    break
                if code.translate(nxt) == code.translate(cur):
                    ps += 1
                else:
                    pn += 1
                cur = nxt
            if ok:
                paths += 1
                paths_s += ps
                paths_n += pn
        if paths:
            sd[i, j] = paths_s / paths
            nd[i, j] = paths_n / paths
            valid[i, j] = True
    return syn_sites, sd, nd, valid


def codon_matrix(sequences: list[str]) -> np.ndarray:
    """``(N, C)`` codon indices (0..63) or ``-1`` for incomplete/ambiguous codons."""

    n_codons = min(len(s) for s in sequences) // 3 if sequences else 0
    out = np.full((len(sequences), n_codons), -1, dtype=np.int64)
    for r, seq in enumerate(sequences):
        s = seq.upper().replace("U", "T")
        for c in range(n_codons):
            out[r, c] = CODON_INDEX.get(s[3 * c : 3 * c + 3], -1)
    return out


def _jc(p: float) -> float | None:
    if p >= 0.75:
        return None
    return -0.75 * math.log(1 - 4 * p / 3)


def pairwise_dnds(
    codons: np.ndarray,
    ids: list[str],
    code_id: int = 1,
    max_pairs: int = DEFAULT_MAX_PAIRS,
    seed: int = 0,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    syn_sites, sd_t, nd_t, valid_t = _tables(code_id)
    n = codons.shape[0]
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
    rows = []
    for i, j in pairs:
        a, b = codons[i], codons[j]
        ok = (a >= 0) & (b >= 0)
        ok[ok] = valid_t[a[ok], b[ok]]
        m = int(ok.sum())
        row: dict[str, Any] = {"sequence_a": ids[i], "sequence_b": ids[j], "codons": m}
        if m < MIN_CODONS:
            rows.append(row | {"dn": None, "ds": None, "dn_ds": None, "note": "too few codons"})
            continue
        s_sites = float((syn_sites[a[ok]] + syn_sites[b[ok]]).sum() / 2)
        n_sites = 3 * m - s_sites
        sd = float(sd_t[a[ok], b[ok]].sum())
        nd = float(nd_t[a[ok], b[ok]].sum())
        ps = sd / s_sites if s_sites else math.nan
        pn = nd / n_sites if n_sites else math.nan
        ds = _jc(ps) if math.isfinite(ps) else None
        dn = _jc(pn) if math.isfinite(pn) else None
        ratio = dn / ds if dn is not None and ds is not None and ds > 0 else None
        note = ""
        if ds is None or dn is None:
            note = "saturated (p >= 0.75)"
        elif ds == 0:
            note = "dS = 0"
        rows.append(
            row
            | {
                "synonymous_sites": s_sites,
                "nonsynonymous_sites": n_sites,
                "synonymous_differences": sd,
                "nonsynonymous_differences": nd,
                "ps": ps,
                "pn": pn,
                "ds": ds,
                "dn": dn,
                "dn_ds": ratio,
                "note": note,
            }
        )
    ratios = [r["dn_ds"] for r in rows if r.get("dn_ds") is not None]
    dss = [r["ds"] for r in rows if r.get("ds") is not None]
    dns = [r["dn"] for r in rows if r.get("dn") is not None]
    summary = {
        "dnds_method": "Nei-Gojobori 1986 with Jukes-Cantor correction",
        "dnds_pair_method": method,
        "dnds_pairs": len(rows),
        "dnds_pairs_defined": len(ratios),
        "median_pairwise_dn_ds": float(np.median(ratios)) if ratios else None,
        "mean_pairwise_ds": float(np.mean(dss)) if dss else None,
        "mean_pairwise_dn": float(np.mean(dns)) if dns else None,
        "dnds_pairs_saturated": sum("saturated" in r.get("note", "") for r in rows),
    }
    return summary, rows


def effective_number_of_codons(codons: np.ndarray, code: GeneticCode) -> float | None:
    counts = np.bincount(codons[codons >= 0].ravel(), minlength=64)
    families: dict[str, list[int]] = {}
    for i, c in enumerate(CODONS):
        aa = code.translate(c)
        if aa != "*":
            families.setdefault(aa, []).append(i)
    by_size: dict[int, list[float]] = {}
    n_by_size: dict[int, int] = {}
    for idx in families.values():
        size = len(idx)
        n_by_size[size] = n_by_size.get(size, 0) + 1
        if size == 1:
            continue
        x = counts[idx].astype(float)
        tot = x.sum()
        if tot < 2:
            continue
        f = (tot * ((x / tot) ** 2).sum() - 1) / (tot - 1)
        by_size.setdefault(size, []).append(f)
    enc = float(n_by_size.get(1, 0))
    means = {s: float(np.mean(v)) for s, v in by_size.items() if v and np.mean(v) > 0}
    for size, n_aa in n_by_size.items():
        if size == 1:
            continue
        f = means.get(size)
        if f is None and size == 3 and 2 in means and 4 in means:
            f = (means[2] + means[4]) / 2  # Wright (1990)
        if f is None:
            return None
        enc += n_aa / f
    return float(min(enc, 61.0))


def codon_evolution_analysis(
    sequences: list[str], ids: list[str], code_id: int = 1, seed: int = 0
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    code = genetic_code(code_id)
    codons = codon_matrix(sequences)
    aa_of = np.array([code.translate(c) for c in CODONS])
    syn_sites = 0
    nonsyn_sites = 0
    variable_codon_sites = 0
    fourfold_pi: list[float] = []
    for c in range(codons.shape[1]):
        col = codons[:, c]
        obs = col[col >= 0]
        obs = obs[aa_of[obs] != "*"]
        if obs.size < 2:
            continue
        uniq = np.unique(obs)
        if uniq.size > 1:
            variable_codon_sites += 1
            if np.unique(aa_of[uniq]).size == 1:
                syn_sites += 1
            else:
                nonsyn_sites += 1
        consensus = CODONS[int(np.bincount(obs, minlength=64).argmax())]
        if all(code.translate(consensus[:2] + b) == code.translate(consensus) for b in BASES):
            thirds = np.array([CODONS[i][2] for i in obs if CODONS[i][:2] == consensus[:2]])
            if thirds.size >= 2:
                _, cnt = np.unique(thirds, return_counts=True)
                m = thirds.size
                fourfold_pi.append(float((m * m - (cnt**2).sum()) / (m * (m - 1))))
    summary: dict[str, Any] = {
        "genetic_code": code.id,
        "codon_sites": int(codons.shape[1]),
        "variable_codon_sites": variable_codon_sites,
        "synonymous_only_codon_sites": syn_sites,
        "nonsynonymous_codon_sites": nonsyn_sites,
        "effective_number_of_codons": effective_number_of_codons(codons, code),
        "fourfold_sites": len(fourfold_pi),
        "fourfold_nucleotide_diversity": float(np.mean(fourfold_pi)) if fourfold_pi else None,
    }
    dnds_summary, rows = pairwise_dnds(codons, ids, code_id, seed=seed)
    summary.update(dnds_summary)
    return summary, {"pairwise_dnds": rows}
