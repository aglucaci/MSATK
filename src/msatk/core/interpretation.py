"""Plain-language interpretation of alignment profiles."""

from __future__ import annotations

from typing import Any


def interpret_summary(summary: dict[str, Any], warnings: list[str]) -> str:
    """Create a concise human-readable interpretation."""

    nseq = int(summary.get("number_of_sequences", 0))
    length = int(summary.get("alignment_length", 0))
    molecule = str(summary.get("molecule_type", "unknown"))
    gap_fraction = float(summary.get("gap_fraction") or 0.0)
    variable = int(summary.get("variable_sites", 0))
    identity = float(summary.get("mean_pairwise_identity") or 0.0)

    completeness = (
        "Overall missingness is low"
        if gap_fraction < 0.1
        else "Overall missingness is moderate"
        if gap_fraction < 0.3
        else "Overall missingness is high"
    )
    suitability = (
        "suggesting the alignment is suitable for most downstream analyses"
        if gap_fraction < 0.1 and identity > 0.7
        else "so flagged sequences and high-gap sites should be reviewed before downstream analysis"
    )
    codon_note = (
        " The alignment appears codon-aware; stop codon and frame-related warnings should be reviewed before selection analyses."
        if molecule == "codon"
        else ""
    )
    warning_note = (
        f" MSATK generated {len([w for w in warnings if w.startswith('WARNING')])} warning(s) that should be inspected."
        if any(w.startswith("WARNING") for w in warnings)
        else " No major QC warnings were triggered by the current thresholds."
    )
    return (
        f"This {molecule} alignment contains {nseq:,} sequences and {length:,} sites, with "
        f"{variable:,} variable sites and mean pairwise identity of {identity:.1%}. "
        f"{completeness} at {gap_fraction:.1%}, {suitability}.{codon_note}{warning_note}"
    )


def recommended_next_steps(summary: dict[str, Any], warnings: list[str]) -> list[str]:
    """Suggest practical next steps from profile results."""

    steps: list[str] = []
    gap_fraction = float(summary.get("gap_fraction") or 0.0)
    molecule = str(summary.get("molecule_type", "unknown"))
    if gap_fraction > 0.1:
        steps.append(
            "Inspect high-gap sequences and consider trimming gap-heavy alignment columns."
        )
    if any(
        "internal stop" in warning.lower() or "stop codon" in warning.lower()
        for warning in warnings
    ):
        steps.append(
            "Review sequences with stop codons for frameshifts, annotation errors, or pseudogenization."
        )
    if molecule == "codon":
        steps.append("Review codon-aware tables before HyPhy, PAML, or other selection analyses.")
    if not steps:
        steps.append(
            "Proceed with downstream phylogenetic, comparative, or modeling analyses as appropriate."
        )
    return steps
