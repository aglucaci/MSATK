"""Unified loading and inspection of every supported alignment type."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from msatk.alignments.blocks import BlockAlignment, parse_axt, parse_maf, parse_xmfa
from msatk.alignments.graphs import GraphAlignment, parse_gfa
from msatk.alignments.mappings import MappingAlignment, parse_blast_tab, parse_paf, parse_psl
from msatk.alignments.profiles import ProfileModel, parse_hhm, parse_hmmer
from msatk.alignments.reads import ReadAlignment, parse_bam_cram, parse_sam
from msatk.alignments.registry import KINDS, get_format
from msatk.exceptions import AlignmentFormatError
from msatk.io.readers import (
    EXTENSION_FORMATS,
    Partition,
    SequenceAnnotations,
    canonical_suffix,
    detect_format,
    is_gzipped,
    parse_partition_file,
    read_sequence_alignment,
    read_text,
)
from msatk.models import Alignment


@dataclass
class SequenceAlignment:
    """A rectangular pairwise or multiple sequence alignment plus its annotations."""

    kind = "sequence"
    alignment: Alignment
    annotations: SequenceAnnotations
    molecule_type: str
    fmt: str
    source: str = ""
    genetic_code: int = 1
    partitions: list[Partition] = field(default_factory=list)

    @property
    def scope(self) -> str:
        n = self.alignment.n_sequences
        return "single" if n == 1 else "pairwise" if n == 2 else "multiple"

    def _run(self, pairwise: str = "auto") -> tuple[Any, ...]:
        from msatk.stats.alignment import compute_statistics
        from msatk.stats.specialized import specialized_analyses

        bundle = compute_statistics(
            self.alignment, self.molecule_type, pairwise=pairwise, genetic_code=self.genetic_code
        )
        return (
            bundle,
            *specialized_analyses(bundle, self.molecule_type, self.annotations, self.partitions),
        )

    def summary(self, pairwise: str = "auto") -> dict[str, Any]:
        bundle, special, _, why = self._run(pairwise)
        return {
            "alignment_kind": "sequence",
            "alignment_scope": self.scope,
            "alignment_method": self.annotations.method,
            **bundle.summary,
            "specialized": special,
            "applicability": why,
        }

    def tables(self, pairwise: str = "auto") -> dict[str, list[dict[str, Any]]]:
        bundle, _, tables, _ = self._run(pairwise)
        return {"sites": bundle.sites, "sequences": bundle.sequences, **tables}


AnyAlignment = (
    SequenceAlignment
    | BlockAlignment
    | ReadAlignment
    | MappingAlignment
    | GraphAlignment
    | ProfileModel
)


def _parse_gaf(text: str, source: str = "") -> MappingAlignment:
    return parse_paf(text, source, fmt="gaf")


_TEXT_PARSERS: dict[str, Callable[[str, str], AnyAlignment]] = {
    "maf": parse_maf,
    "xmfa": parse_xmfa,
    "axt": parse_axt,
    "sam": parse_sam,
    "paf": parse_paf,
    "gaf": _parse_gaf,
    "psl": parse_psl,
    "blast_tab": parse_blast_tab,
    "gfa": parse_gfa,
    "hmmer": parse_hmmer,
    "hhm": parse_hhm,
}


def detect_input_format(path: str | Path) -> str:
    source = Path(path)
    suffix = canonical_suffix(source)
    if suffix in {"bam", "cram"}:
        return suffix
    by_ext = EXTENSION_FORMATS.get(suffix)
    if by_ext in {"gbz", "odgi", "hal", "pdb", "mmcif", "chain"}:
        return by_ext
    head = read_text(source)[:200_000]
    return detect_format(head, source)


def load_alignment(
    path: str | Path,
    fmt: str = "auto",
    molecule_type: str = "auto",
    genetic_code: int = 1,
    partitions: str | Path | None = None,
    reference: str | None = None,
    validation_mode: str = "permissive",
) -> AnyAlignment:
    """Load any supported alignment and return the matching data model.

    The returned object has ``.kind`` (``sequence``, ``block``, ``read``, ``mapping``,
    ``graph`` or ``profile_model``), ``.summary()`` and ``.tables()``.
    """

    source = Path(path)
    name = detect_input_format(source) if fmt == "auto" else EXTENSION_FORMATS.get(fmt, fmt)
    cap = get_format(name)
    if cap.level == "planned":
        raise AlignmentFormatError(
            f"{cap.label} is recognised but not supported yet ({cap.notes or 'planned'})."
        )
    if cap.kind == "sequence":
        from msatk.core.detect import normalize_requested_type

        alignment, ann = read_sequence_alignment(source, name, validation_mode)
        parts = parse_partition_file(partitions, alignment.length) if partitions else []
        return SequenceAlignment(
            alignment,
            ann,
            normalize_requested_type(molecule_type, alignment),
            name,
            str(source),
            genetic_code,
            parts,
        )
    if name in {"bam", "cram"}:
        return parse_bam_cram(source, name, reference)
    parser = _TEXT_PARSERS.get(name)
    if parser is None:
        raise AlignmentFormatError(f"No reader for format '{name}'.")
    return parser(read_text(source), str(source))


def inspect_alignment(path: str | Path, fmt: str = "auto", **kwargs: Any) -> dict[str, Any]:
    """Describe what an input contains and which analyses are valid for it."""

    source = Path(path)
    name = detect_input_format(source) if fmt == "auto" else EXTENSION_FORMATS.get(fmt, fmt)
    cap = get_format(name)
    info: dict[str, Any] = {
        "path": str(source),
        "compressed": "gzip" if is_gzipped(source) else None,
        "format": name,
        "format_label": cap.label,
        "format_family": cap.family,
        "support_level": cap.level,
        "kind": cap.kind,
        "data_model": KINDS[cap.kind]["model"],
        "preserved_annotations": list(cap.preserves),
        "valid_analyses": KINDS[cap.kind]["analyses"],
        "not_computed_by_design": KINDS[cap.kind]["not_valid"],
    }
    if cap.level == "planned":
        info["status"] = f"not readable yet: {cap.notes or 'planned'}"
        return info
    obj = load_alignment(source, name, **kwargs)
    if isinstance(obj, SequenceAlignment):
        from msatk.stats.alignment import compute_statistics
        from msatk.stats.specialized import specialized_analyses

        bundle = compute_statistics(
            obj.alignment, obj.molecule_type, pairwise="none", genetic_code=obj.genetic_code
        )
        _, _, why = specialized_analyses(bundle, obj.molecule_type, obj.annotations, obj.partitions)
        info.update(
            {
                "molecule_type": obj.molecule_type,
                "scope": obj.scope,
                "method": obj.annotations.method,
                "sequences": obj.alignment.n_sequences,
                "columns": obj.alignment.length,
                "rectangular": obj.alignment.is_rectangular,
                "annotations_found": obj.annotations.present()
                + (["partitions (file)"] if obj.partitions else []),
                "specialized_analyses": why,
            }
        )
    else:
        summary = obj.summary()
        info["contents"] = {
            k: v for k, v in summary.items() if k not in {"alignment_kind", "detected_format"}
        }
    info["status"] = "readable"
    return info
