"""Reference-based read alignments: SAM (native), BAM/CRAM (via pysam).

Reads are summarised with mapping, CIGAR, quality and coverage statistics, not
as rows of a multiple sequence alignment: an all-vs-all identity matrix of mapped
reads is neither meaningful nor tractable.

Depth counts aligned bases (CIGAR ``M``, ``=``, ``X``) of primary, non-duplicate,
mapped reads, matching the default of ``samtools depth``. Coverage is computed
only for references with a declared length (``@SQ LN``) and when the total
reference length is at most ``max_coverage_bases``.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean, median
from typing import Any

import numpy as np

from msatk.exceptions import AlignmentFormatError, MissingDependencyError

FLAG_PAIRED, FLAG_PROPER, FLAG_UNMAPPED = 0x1, 0x2, 0x4
FLAG_READ1, FLAG_SECONDARY, FLAG_QCFAIL, FLAG_DUP, FLAG_SUPPLEMENTARY = (
    0x40,
    0x100,
    0x200,
    0x400,
    0x800,
)
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_OPS = "MIDNSHP=X"


@dataclass
class ReadRecord:
    name: str
    flag: int
    reference: str | None
    pos: int  # 0-based leftmost aligned position (-1 if unmapped)
    mapq: int
    cigar: list[tuple[str, int]]
    length: int
    tlen: int = 0
    nm: int | None = None
    mean_quality: float | None = None

    @property
    def mapped(self) -> bool:
        return not self.flag & FLAG_UNMAPPED and self.reference is not None

    @property
    def primary(self) -> bool:
        return not self.flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)


@dataclass
class ReadAlignment:
    kind = "read"
    records: list[ReadRecord]
    references: dict[str, int]
    fmt: str
    source: str = ""
    meta: dict[str, Any] = field(default_factory=dict)

    def _cigar_totals(self, reads: list[ReadRecord]) -> Counter[str]:
        totals: Counter[str] = Counter()
        for r in reads:
            for op, n in r.cigar:
                totals[op] += n
        return totals

    def coverage(self, max_coverage_bases: float = 5e8) -> dict[str, np.ndarray] | None:
        if not self.references or sum(self.references.values()) > max_coverage_bases:
            return None
        diff = {name: np.zeros(n + 1, dtype=np.int32) for name, n in self.references.items()}
        for r in self.records:
            if not (r.mapped and r.primary) or r.flag & FLAG_DUP or r.reference not in diff:
                continue
            d = diff[r.reference]
            pos = r.pos
            for op, n in r.cigar:
                if op in "M=X":
                    d[min(pos, len(d) - 1)] += 1
                    d[min(pos + n, len(d) - 1)] -= 1
                    pos += n
                elif op in "DN":
                    pos += n
        return {name: np.cumsum(d)[:-1] for name, d in diff.items()}

    def reference_rows(self) -> list[dict[str, Any]]:
        depth = self.coverage()
        counts = Counter(r.reference for r in self.records if r.mapped and r.primary)
        rows = []
        for name, length in self.references.items():
            row: dict[str, Any] = {
                "reference": name,
                "length": length,
                "primary_mapped_reads": counts.get(name, 0),
            }
            if depth is not None:
                dpt = depth[name]
                row.update(
                    {
                        "mean_depth": float(dpt.mean()) if len(dpt) else None,
                        "median_depth": float(np.median(dpt)) if len(dpt) else None,
                        "breadth_1x": float((dpt >= 1).mean()) if len(dpt) else None,
                        "breadth_10x": float((dpt >= 10).mean()) if len(dpt) else None,
                        "breadth_30x": float((dpt >= 30).mean()) if len(dpt) else None,
                    }
                )
            rows.append(row)
        return rows

    def summary(self) -> dict[str, Any]:
        recs = self.records
        primary = [r for r in recs if r.primary]
        mapped = [r for r in primary if r.mapped]
        mapq = [r.mapq for r in mapped if r.mapq != 255]
        cig = self._cigar_totals(mapped)
        aligned = cig["M"] + cig["="] + cig["X"]
        nm_known = [r for r in mapped if r.nm is not None]
        mismatch = None
        if cig["="] + cig["X"] and not cig["M"]:
            mismatch = cig["X"] / aligned if aligned else None
        elif mapped and len(nm_known) == len(mapped):
            mism = sum(r.nm or 0 for r in mapped) - cig["I"] - cig["D"]
            mismatch = max(mism, 0) / aligned if aligned else None
        quals = [r.mean_quality for r in primary if r.mean_quality is not None]
        lengths = sorted(r.length for r in primary if r.length)
        tlens = [
            abs(r.tlen) for r in mapped if r.flag & FLAG_PROPER and r.flag & FLAG_READ1 and r.tlen
        ]
        refs = self.reference_rows()
        total_len = sum(self.references.values())
        cov = [r for r in refs if "mean_depth" in r]
        return {
            "alignment_kind": "read",
            "detected_format": self.fmt,
            "references": len(self.references),
            "reference_bases": total_len or None,
            "records": len(recs),
            "primary_reads": len(primary),
            "secondary_records": sum(1 for r in recs if r.flag & FLAG_SECONDARY),
            "supplementary_records": sum(1 for r in recs if r.flag & FLAG_SUPPLEMENTARY),
            "unmapped_reads": len(primary) - len(mapped),
            "mapping_rate": len(mapped) / len(primary) if primary else None,
            "duplicate_reads": sum(1 for r in primary if r.flag & FLAG_DUP),
            "qc_fail_reads": sum(1 for r in primary if r.flag & FLAG_QCFAIL),
            "paired_reads": sum(1 for r in primary if r.flag & FLAG_PAIRED),
            "properly_paired_reads": sum(1 for r in primary if r.flag & FLAG_PROPER),
            "mean_mapq": mean(mapq) if mapq else None,
            "median_mapq": median(mapq) if mapq else None,
            "mapq_ge_30_fraction": sum(q >= 30 for q in mapq) / len(mapq) if mapq else None,
            "mapq_0_fraction": sum(q == 0 for q in mapq) / len(mapq) if mapq else None,
            "aligned_bases": aligned,
            "inserted_bases": cig["I"],
            "deleted_bases": cig["D"],
            "indel_rate": (cig["I"] + cig["D"]) / aligned if aligned else None,
            "mismatch_rate": mismatch,
            "soft_clipped_bases": cig["S"],
            "hard_clipped_bases": cig["H"],
            "clipped_read_fraction": sum(1 for r in mapped if any(op in "SH" for op, _ in r.cigar))
            / len(mapped)
            if mapped
            else None,
            "spliced_reads": sum(1 for r in mapped if any(op == "N" for op, _ in r.cigar)),
            "mean_read_length": mean(lengths) if lengths else None,
            "n50_read_length": _n50(lengths),
            "mean_base_quality": mean(quals) if quals else None,
            "median_insert_size": median(tlens) if tlens else None,
            "mean_depth": sum(r["mean_depth"] * r["length"] for r in cov) / total_len
            if cov and total_len
            else None,
            "breadth_1x": sum(r["breadth_1x"] * r["length"] for r in cov) / total_len
            if cov and total_len
            else None,
            "breadth_10x": sum(r["breadth_10x"] * r["length"] for r in cov) / total_len
            if cov and total_len
            else None,
        }

    def tables(self) -> dict[str, list[dict[str, Any]]]:
        hist = Counter(r.mapq for r in self.records if r.mapped and r.primary)
        return {
            "references": self.reference_rows(),
            "mapq_histogram": [{"mapq": q, "reads": n} for q, n in sorted(hist.items())],
        }


def _n50(values: list[int]) -> int | None:
    if not values:
        return None
    half = sum(values) / 2
    acc = 0
    for v in sorted(values, reverse=True):
        acc += v
        if acc >= half:
            return v
    return None


def parse_sam(text: str, source: str = "") -> ReadAlignment:
    refs: dict[str, int] = {}
    records: list[ReadRecord] = []
    for raw in text.splitlines():
        if not raw.strip():
            continue
        if raw.startswith("@"):
            if raw.startswith("@SQ"):
                tags = dict(t.split(":", 1) for t in raw.split("\t")[1:] if ":" in t)
                if "SN" in tags and "LN" in tags:
                    refs[tags["SN"]] = int(tags["LN"])
            continue
        p = raw.rstrip("\n").split("\t")
        if len(p) < 11:
            raise AlignmentFormatError(f"Invalid SAM record: {raw[:80]}")
        flag = int(p[1])
        cigar = [(op, int(n)) for n, op in _CIGAR.findall(p[5])] if p[5] != "*" else []
        nm = next((int(t[5:]) for t in p[11:] if t.startswith("NM:i:")), None)
        qual = (
            None if p[10] == "*" else float(np.mean(np.frombuffer(p[10].encode(), np.uint8) - 33))
        )
        seq_len = len(p[9]) if p[9] != "*" else sum(n for op, n in cigar if op in "MIS=X")
        records.append(
            ReadRecord(
                name=p[0],
                flag=flag,
                reference=None if p[2] == "*" else p[2],
                pos=int(p[3]) - 1,
                mapq=int(p[4]),
                cigar=cigar,
                length=seq_len,
                tlen=int(p[8]),
                nm=nm,
                mean_quality=qual,
            )
        )
    if not records:
        raise AlignmentFormatError("No SAM alignment records found.")
    return ReadAlignment(records, refs, "sam", source)


def parse_bam_cram(path: str | Path, fmt: str, reference: str | None = None) -> ReadAlignment:
    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise MissingDependencyError(
            f"{fmt.upper()} input requires pysam: `pip install pysam` or "
            "`mamba install -c bioconda pysam`."
        ) from exc
    mode = "rc" if fmt == "cram" else "rb"
    kwargs = {"reference_filename": reference} if reference else {}
    records: list[ReadRecord] = []
    with pysam.AlignmentFile(str(path), mode, **kwargs) as handle:
        refs = dict(zip(handle.references, handle.lengths))
        for r in handle.fetch(until_eof=True):
            quals = r.query_qualities
            records.append(
                ReadRecord(
                    name=r.query_name,
                    flag=r.flag,
                    reference=None if r.is_unmapped else r.reference_name,
                    pos=r.reference_start if not r.is_unmapped else -1,
                    mapq=r.mapping_quality,
                    cigar=[(_OPS[op], n) for op, n in (r.cigartuples or [])],
                    length=r.query_length or r.infer_read_length() or 0,
                    tlen=r.template_length,
                    nm=r.get_tag("NM") if r.has_tag("NM") else None,
                    mean_quality=float(np.mean(quals))
                    if quals is not None and len(quals)
                    else None,
                )
            )
    return ReadAlignment(records, refs, fmt, str(path))
