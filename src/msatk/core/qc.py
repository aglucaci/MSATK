"""Quality-control rules producing severity-ranked, explainable findings."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from msatk.core.genetic_code import genetic_code as _genetic_code
from msatk.models import Alignment

SEVERITY_ORDER = {"error": 0, "warning": 1, "info": 2, "ok": 3}


@dataclass
class QCFinding:
    code: str
    severity: str
    scope: str
    message: str
    count: int = 0
    suggestion: str = ""
    items: list[str] = field(default_factory=list)

    def as_text(self) -> str:
        prefix = {"error": "ERROR", "warning": "WARNING", "info": "NOTE", "ok": "OK"}[self.severity]
        return f"{prefix}: {self.message}"

    def as_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["items"] = ";".join(self.items[:200])
        return data


@dataclass
class QCThresholds:
    max_gap_seq: float = 0.3
    max_gap_site: float = 0.5
    min_consensus_identity: float = 0.5
    robust_z: float = 3.5
    max_ambiguous_seq: float = 0.05


def flagged_sequences(
    per_sequence: list[dict[str, Any]], max_gap_seq: float = 0.3
) -> list[dict[str, Any]]:
    return [row for row in per_sequence if float(row.get("gap_fraction") or 0.0) > max_gap_seq]


def flagged_sites(
    per_site: list[dict[str, Any]], max_gap_site: float = 0.5
) -> list[dict[str, Any]]:
    return [row for row in per_site if float(row.get("gap_fraction") or 0.0) > max_gap_site]


def qc_findings(
    alignment: Alignment,
    molecule_type: str,
    per_sequence: list[dict[str, Any]],
    per_site: list[dict[str, Any]],
    summary: dict[str, Any] | None = None,
    thresholds: QCThresholds | None = None,
    genetic_code: int | None = None,
) -> list[QCFinding]:
    t = thresholds or QCThresholds()
    summary = summary or {}
    code = genetic_code or int(summary.get("genetic_code") or 1)
    stops = _genetic_code(code).stops
    out: list[QCFinding] = []

    if not alignment.is_rectangular:
        out.append(
            QCFinding(
                "UNEQUAL_LENGTHS",
                "warning",
                "alignment",
                f"Sequence lengths are inconsistent; {alignment.padded_cells} cell(s) were padded "
                "as missing data ('?') to form a rectangle. Padding is not treated as gaps.",
                alignment.padded_cells,
                "Confirm the input is truly aligned; realign or use --strict to reject it.",
            )
        )
    dup_ids = sorted({i for i in alignment.ids if alignment.ids.count(i) > 1})
    if dup_ids:
        out.append(
            QCFinding(
                "DUPLICATE_IDS",
                "warning",
                "sequence",
                f"{len(dup_ids)} sequence identifier(s) occur more than once.",
                len(dup_ids),
                "Rename duplicates; many downstream tools require unique identifiers.",
                dup_ids,
            )
        )
    dups = [r["sequence_id"] for r in per_sequence if r.get("duplicate_of")]
    if dups:
        out.append(
            QCFinding(
                "DUPLICATE_SEQUENCES",
                "info",
                "sequence",
                f"{len(dups)} sequence(s) are identical to an earlier sequence.",
                len(dups),
                "Consider deduplicating before phylogenetic or ML analyses to avoid redundancy bias.",
                dups,
            )
        )
    if summary.get("invalid_cells"):
        out.append(
            QCFinding(
                "INVALID_CHARACTERS",
                "warning",
                "alignment",
                f"{summary['invalid_cells']} cell(s) contain characters outside the "
                f"{molecule_type} alphabet.",
                int(summary["invalid_cells"]),
                "Check the molecule type (--type) and input encoding.",
            )
        )
    seq_flags = flagged_sequences(per_sequence, t.max_gap_seq)
    if seq_flags:
        out.append(
            QCFinding(
                "HIGH_GAP_SEQUENCES",
                "warning",
                "sequence",
                f"{len(seq_flags)} sequences have >{t.max_gap_seq:.0%} gaps.",
                len(seq_flags),
                "Inspect partial or fragmentary sequences; consider removing them.",
                [str(r["sequence_id"]) for r in seq_flags],
            )
        )
    site_flags = flagged_sites(per_site, t.max_gap_site)
    if site_flags:
        out.append(
            QCFinding(
                "HIGH_GAP_SITES",
                "warning",
                "site",
                f"{len(site_flags)} alignment columns have >{t.max_gap_site:.0%} gaps.",
                len(site_flags),
                "Consider trimming gap-rich columns (e.g. trimAl, ClipKIT) before tree inference.",
                [str(r["site_index"]) for r in site_flags],
            )
        )
    unobserved = [r for r in per_site if r.get("site_class") == "unobserved"]
    if unobserved:
        out.append(
            QCFinding(
                "UNOBSERVED_SITES",
                "info",
                "site",
                f"{len(unobserved)} column(s) contain no canonical residue; their conservation "
                "and entropy are undefined and they are not counted as constant sites.",
                len(unobserved),
                "Remove empty columns before analysis.",
                [str(r["site_index"]) for r in unobserved],
            )
        )
    divergent = [
        r
        for r in per_sequence
        if (
            r.get("consensus_identity") is not None
            and r["consensus_identity"] < t.min_consensus_identity
        )
        or (
            r.get("consensus_identity_robust_z") is not None
            and r["consensus_identity_robust_z"] < -t.robust_z
        )
    ]
    if divergent:
        out.append(
            QCFinding(
                "DIVERGENT_SEQUENCES",
                "warning",
                "sequence",
                f"{len(divergent)} sequences have unusually low identity to the consensus.",
                len(divergent),
                "Check for misalignment, contamination, paralogy or reverse-complemented input.",
                [str(r["sequence_id"]) for r in divergent],
            )
        )
    comp = [
        r
        for r in per_sequence
        if r.get("composition_robust_z") is not None and r["composition_robust_z"] > t.robust_z
    ]
    if comp:
        out.append(
            QCFinding(
                "COMPOSITION_OUTLIERS",
                "warning",
                "sequence",
                f"{len(comp)} sequences have atypical residue composition.",
                len(comp),
                "Compositional heterogeneity can bias phylogenetic inference; inspect these rows.",
                [str(r["sequence_id"]) for r in comp],
            )
        )
    amb = [r for r in per_sequence if float(r.get("ambiguous_fraction") or 0) > t.max_ambiguous_seq]
    if amb:
        out.append(
            QCFinding(
                "HIGH_AMBIGUITY_SEQUENCES",
                "warning",
                "sequence",
                f"{len(amb)} sequences have >{t.max_ambiguous_seq:.0%} ambiguity codes.",
                len(amb),
                "Low-quality or mixed samples; consider masking or removal.",
                [str(r["sequence_id"]) for r in amb],
            )
        )
    no_comp = int(summary.get("pairwise_pairs_without_comparable_sites") or 0)
    if no_comp:
        out.append(
            QCFinding(
                "NON_OVERLAPPING_PAIRS",
                "warning",
                "pair",
                f"{no_comp} sequence pair(s) share no comparable site; their identity is undefined.",
                no_comp,
                "Non-overlapping fragments cannot be compared directly; distance-based methods may fail.",
            )
        )
    if molecule_type == "codon":
        out.extend(_codon_findings(alignment, stops))
    if molecule_type in {"protein", "translated_cds"}:
        stop_rows = [r.id for r in alignment.records if "*" in r.sequence.rstrip("-.*")]
        if stop_rows:
            out.append(
                QCFinding(
                    "INTERNAL_STOPS_PROTEIN",
                    "warning",
                    "sequence",
                    f"{len(stop_rows)} protein sequence(s) contain internal stop symbols ('*').",
                    len(stop_rows),
                    "Check translation frame and genetic code.",
                    stop_rows,
                )
            )
    if alignment.n_sequences == 1:
        out.append(
            QCFinding(
                "SINGLE_SEQUENCE",
                "info",
                "alignment",
                "Only one sequence: variability, pairwise and diversity statistics are undefined.",
                1,
            )
        )
    if not [f for f in out if f.severity in {"error", "warning"}]:
        out.append(
            QCFinding(
                "NO_MAJOR_ISSUES",
                "ok",
                "alignment",
                "No major MSATK QC warnings triggered by current thresholds.",
            )
        )
    out.sort(key=lambda f: (SEVERITY_ORDER[f.severity], -f.count))
    return out


def _codon_findings(alignment: Alignment, stops: frozenset[str]) -> list[QCFinding]:
    out: list[QCFinding] = []
    if alignment.length % 3 != 0:
        out.append(
            QCFinding(
                "CODON_LENGTH_NOT_MULTIPLE_OF_3",
                "warning",
                "alignment",
                "Alignment length is not divisible by 3; codon-aware analysis may be invalid.",
                1,
                "Use a codon-aware aligner (e.g. MACSE, PAL2NAL) or analyse as --type dna.",
            )
        )
    stop_ids: list[str] = []
    partial = 0
    for record in alignment.records:
        seq = record.sequence.upper().replace("U", "T")
        codons = [seq[i : i + 3] for i in range(0, len(seq) - 2, 3)]
        internal = [c for c in codons[:-1] if c in stops]
        if internal:
            stop_ids.append(record.id)
            out.append(
                QCFinding(
                    "INTERNAL_STOP_CODON",
                    "warning",
                    "sequence",
                    f"Sequence {record.id} contains {len(internal)} internal stop codon(s).",
                    len(internal),
                    "Check for frameshifts, sequencing errors or pseudogenes before selection analyses.",
                    [record.id],
                )
            )
        partial += sum(1 for c in codons if 0 < sum(ch in "-." for ch in c) < 3)
    if partial:
        out.append(
            QCFinding(
                "PARTIAL_CODON_GAPS",
                "warning",
                "site",
                f"{partial} codon(s) are partially gapped, indicating frame-disrupting indels.",
                partial,
                "Codon models assume gaps in multiples of three within codons.",
            )
        )
    return out


def qc_warnings(
    alignment: Alignment,
    molecule_type: str,
    per_sequence: list[dict[str, Any]],
    per_site: list[dict[str, Any]],
    max_gap_seq: float = 0.3,
    max_gap_site: float = 0.5,
    summary: dict[str, Any] | None = None,
) -> list[str]:
    """Text form of :func:`qc_findings` (kept for backward compatibility)."""

    findings = qc_findings(
        alignment,
        molecule_type,
        per_sequence,
        per_site,
        summary=summary,
        thresholds=QCThresholds(max_gap_seq=max_gap_seq, max_gap_site=max_gap_site),
    )
    return [f.as_text() for f in findings]
