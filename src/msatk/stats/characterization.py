"""Alignment characterization: named analyses grouped into metric profiles.

Each analysis answers one question and runs only when it is valid for the input;
``applicability`` records why it ran or was skipped. Profiles select analyses:

==============  =========================================================================
``core``        completeness, gap structure, redundancy (N_eff)
``comprehensive`` core + site patterns, composition heterogeneity, windows,
                protein property variation (proteins)
``evolutionary`` comprehensive + symmetry tests, codon evolution (codon alignments)
``coevolution`` comprehensive + mutual information with APC
``all``         every analysis
==============  =========================================================================

Divergence from the consensus and its robust z-score, nucleotide diversity,
Watterson's theta, Tajima's D and haplotype diversity are part of the base
statistics (``per_sequence_stats``, ``summary``) and are always computed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from msatk.stats.alignment import StatisticsBundle
from msatk.stats.pairwise import PairwiseResult, exact_pairwise

ANALYSES = (
    "completeness",
    "gap_structure",
    "redundancy",
    "site_patterns",
    "composition_heterogeneity",
    "windows",
    "protein_properties",
    "symmetry_tests",
    "codon_evolution",
    "coevolution",
)
_CORE = ("completeness", "gap_structure", "redundancy")
_COMPREHENSIVE = _CORE + (
    "site_patterns",
    "composition_heterogeneity",
    "windows",
    "protein_properties",
)
PROFILES: dict[str, tuple[str, ...]] = {
    "none": (),
    "core": _CORE,
    "comprehensive": _COMPREHENSIVE,
    "evolutionary": _COMPREHENSIVE + ("symmetry_tests", "codon_evolution"),
    "coevolution": _COMPREHENSIVE + ("coevolution",),
    "all": ANALYSES,
}
TITLES = {
    "completeness": "Completeness and pairwise overlap",
    "gap_structure": "Gap structure",
    "redundancy": "Redundancy and effective sequences",
    "site_patterns": "Site diversity and site patterns",
    "composition_heterogeneity": "Composition heterogeneity",
    "windows": "Windowed statistics",
    "protein_properties": "Protein property variation",
    "symmetry_tests": "Matched-pairs symmetry tests",
    "codon_evolution": "Codon evolution",
    "coevolution": "Coevolution (mutual information)",
}
PROTEIN = {"protein", "translated_cds"}


@dataclass
class Characterization:
    profile: str
    summaries: dict[str, dict[str, Any]] = field(default_factory=dict)
    tables: dict[str, list[dict[str, Any]]] = field(default_factory=dict)
    applicability: dict[str, str] = field(default_factory=dict)
    matrices: dict[str, np.ndarray] = field(default_factory=dict)


def resolve(profile: str | None, analyses: list[str] | tuple[str, ...] | None) -> tuple[str, ...]:
    if analyses:
        unknown = [a for a in analyses if a not in ANALYSES]
        if unknown:
            raise ValueError(f"Unknown analyses {unknown}; choose from {ANALYSES}.")
        return tuple(a for a in ANALYSES if a in analyses)
    name = profile or "core"
    if name not in PROFILES:
        raise ValueError(f"Unknown metrics profile {name!r}; choose from {tuple(PROFILES)}.")
    return PROFILES[name]


def characterize(
    bundle: StatisticsBundle,
    molecule_type: str,
    profile: str | None = "core",
    analyses: list[str] | tuple[str, ...] | None = None,
    max_exact_sequences: int = 1000,
    genetic_code: int = 1,
    seed: int = 0,
    options: dict[str, dict[str, Any]] | None = None,
) -> Characterization:
    """Run the analyses of ``profile`` (or the explicit ``analyses``) on ``bundle``."""

    selected = resolve(profile, analyses)
    opts = options or {}
    enc = bundle.encoded
    n = enc.n
    result = Characterization(profile="custom" if analyses else (profile or "core"))
    why = result.applicability
    supported = molecule_type not in {"mixed", "unknown"}

    pw: PairwiseResult | None = bundle.pairwise if bundle.pairwise is not None else None
    if pw is not None and not pw.is_exact:
        pw = None
    need_exact = {"redundancy", "site_patterns", "coevolution"} & set(selected)
    if pw is None and need_exact and 2 <= n <= max_exact_sequences:
        pw = exact_pairwise(enc)
    weights = None

    for name in ANALYSES:
        if name not in selected:
            why[name] = "not requested (profile " + result.profile + ")"
            continue
        if not supported:
            why[name] = f"skipped: molecule type {molecule_type!r} has no defined alphabet"
            continue
        o = dict(opts.get(name, {}))
        if name == "completeness":
            from msatk.stats.completeness import completeness_analysis

            s, t, m = completeness_analysis(enc, exact_pairs=2 <= n <= max_exact_sequences, **o)
            if m is not None:
                result.matrices["pair_completeness"] = m
            why[name] = "run" + (
                ""
                if m is not None
                else f" (pairwise overlap skipped: N > max_exact_sequences={max_exact_sequences})"
                if n >= 2
                else ""
            )
        elif name == "gap_structure":
            from msatk.stats.completeness import gap_structure_analysis

            s, t = gap_structure_analysis(enc)
            why[name] = "run"
        elif name == "redundancy":
            if n < 2:
                why[name] = "skipped: fewer than 2 sequences"
                continue
            if pw is None:
                why[name] = (
                    f"skipped: N = {n} exceeds max_exact_sequences = {max_exact_sequences} "
                    "(raise it to compute N_eff)"
                )
                continue
            from msatk.stats.redundancy import redundancy_analysis

            s, t, weights = redundancy_analysis(enc, pw, **o)
            why[name] = f"run: identity threshold {s['identity_threshold']}"
        elif name == "site_patterns":
            from msatk.stats.site_patterns import site_pattern_analysis

            if weights is None and pw is not None:
                from msatk.stats.redundancy import sequence_weights

                weights = sequence_weights(enc, pw)
            s, t = site_pattern_analysis(enc, weights)
            why[name] = "run" + ("" if weights is not None else " (unweighted: no exact pairwise)")
        elif name == "composition_heterogeneity":
            if n < 2:
                why[name] = "skipped: fewer than 2 sequences"
                continue
            from msatk.stats.heterogeneity import composition_heterogeneity

            s, t = composition_heterogeneity(enc)
            why[name] = "run"
        elif name == "windows":
            if enc.width < 20:
                why[name] = f"skipped: alignment has {enc.width} sites (needs >= 20)"
                continue
            from msatk.stats.windows import windowed_analysis

            s, t = windowed_analysis(enc, **o)
            why[name] = f"run: window {s['window']}, step {s['step']}"
        elif name == "protein_properties":
            if molecule_type not in PROTEIN:
                why[name] = "skipped: applies to protein alignments"
                continue
            from msatk.protein.properties import protein_property_analysis

            s, t = protein_property_analysis(enc)
            why[name] = "run: protein alignment"
        elif name == "symmetry_tests":
            if n < 2:
                why[name] = "skipped: fewer than 2 sequences"
                continue
            from msatk.stats.heterogeneity import symmetry_tests

            s, t = symmetry_tests(enc, seed=seed, **o)
            why[name] = f"run: {s['pairs_tested']} pairs ({s['pair_method']})"
        elif name == "codon_evolution":
            if molecule_type != "codon":
                why[name] = "skipped: applies to codon (coding nucleotide) alignments"
                continue
            if n < 2:
                why[name] = "skipped: fewer than 2 sequences"
                continue
            from msatk.codon.evolution import codon_evolution_analysis

            seqs = [r.sequence for r in enc.alignment.records]
            s, t = codon_evolution_analysis(seqs, enc.alignment.ids, genetic_code, seed=seed)
            why[name] = f"run: genetic code {genetic_code}"
        elif name == "coevolution":
            min_n = int(o.pop("min_sequences", 20))
            if n < min_n:
                why[name] = f"skipped: {n} sequences (mutual information needs >= {min_n})"
                continue
            if weights is None and pw is not None:
                from msatk.stats.redundancy import sequence_weights

                weights = sequence_weights(enc, pw)
            from msatk.stats.coevolution import coevolution_analysis

            s, t = coevolution_analysis(enc, weights, **o)
            why[name] = f"run: {s['candidate_sites']} candidate columns" + (
                "" if weights is not None else " (unweighted: no exact pairwise)"
            )
        else:  # pragma: no cover
            continue
        result.summaries[name] = s
        result.tables.update(t)
    return result
