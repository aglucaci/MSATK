"""Pairwise mapping and hit tables: PAF, PSL, BLAST tabular (outfmt 6) and GAF.

Each record is one local alignment between a query and a target (a sequence, or a
graph path for GAF). Summaries cover hit counts, identity, alignment length,
mapping quality, multi-mapping and query/target coverage (union of aligned
intervals over the declared sequence length, when the format provides it).

Identity definitions follow each format: PAF/GAF ``matches / block length``; PSL
``(matches + repMatches) / (matches + repMatches + misMatches)``; BLAST
``pident / 100``.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from statistics import mean, median
from typing import Any

from msatk.exceptions import AlignmentFormatError


@dataclass
class MappingRecord:
    query: str
    query_length: int | None
    query_start: int  # 0-based, half-open
    query_end: int
    strand: str
    target: str
    target_length: int | None
    target_start: int
    target_end: int
    identity: float | None
    alignment_length: int
    mapq: int | None = None
    evalue: float | None = None
    bitscore: float | None = None
    gap_openings: int | None = None


@dataclass
class MappingAlignment:
    kind = "mapping"
    records: list[MappingRecord]
    fmt: str
    source: str = ""
    meta: dict[str, Any] = field(default_factory=dict)

    def _coverage(self, side: str) -> list[dict[str, Any]]:
        spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
        lengths: dict[str, int | None] = {}
        hits: dict[str, int] = defaultdict(int)
        for r in self.records:
            name = getattr(r, side)
            s, e = getattr(r, f"{side}_start"), getattr(r, f"{side}_end")
            spans[name].append((min(s, e), max(s, e)))
            lengths[name] = getattr(r, f"{side}_length")
            hits[name] += 1
        rows = []
        for name, iv in spans.items():
            covered = _union(iv)
            n = lengths.get(name)
            rows.append(
                {
                    side: name,
                    "length": n,
                    "hits": hits[name],
                    "covered_bases": covered,
                    "coverage": covered / n if n else None,
                }
            )
        return rows

    def summary(self) -> dict[str, Any]:
        recs = self.records
        ident = [r.identity for r in recs if r.identity is not None]
        alen = [r.alignment_length for r in recs]
        weighted = (
            sum((r.identity or 0) * r.alignment_length for r in recs if r.identity is not None)
            / sum(r.alignment_length for r in recs if r.identity is not None)
            if ident
            else None
        )
        mapq = [r.mapq for r in recs if r.mapq is not None and r.mapq != 255]
        per_query: dict[str, int] = defaultdict(int)
        for r in recs:
            per_query[r.query] += 1
        q_cov = [row["coverage"] for row in self._coverage("query") if row["coverage"] is not None]
        t_cov = [row["coverage"] for row in self._coverage("target") if row["coverage"] is not None]
        ev = [r.evalue for r in recs if r.evalue is not None]
        return {
            "alignment_kind": "mapping",
            "detected_format": self.fmt,
            "records": len(recs),
            "queries": len(per_query),
            "targets": len({r.target for r in recs}),
            "multi_hit_queries": sum(1 for v in per_query.values() if v > 1),
            "mean_identity": mean(ident) if ident else None,
            "median_identity": median(ident) if ident else None,
            "length_weighted_identity": weighted,
            "total_aligned_length": sum(alen),
            "median_alignment_length": median(alen) if alen else None,
            "minus_strand_fraction": sum(r.strand == "-" for r in recs) / len(recs)
            if recs
            else None,
            "mean_mapq": mean(mapq) if mapq else None,
            "mapq_ge_30_fraction": sum(q >= 30 for q in mapq) / len(mapq) if mapq else None,
            "mean_query_coverage": mean(q_cov) if q_cov else None,
            "mean_target_coverage": mean(t_cov) if t_cov else None,
            "min_evalue": min(ev) if ev else None,
        }

    def tables(self) -> dict[str, list[dict[str, Any]]]:
        return {"queries": self._coverage("query"), "targets": self._coverage("target")}


def _union(iv: list[tuple[int, int]]) -> int:
    iv = sorted(iv)
    total = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s > ce:
            total += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return total + ce - cs


def _tag(fields: list[str], key: str) -> str | None:
    return next((f.split(":", 2)[2] for f in fields if f.startswith(key + ":")), None)


def parse_paf(text: str, source: str = "", fmt: str = "paf") -> MappingAlignment:
    """PAF; with ``fmt='gaf'`` the target is a graph path (column 6)."""

    out = []
    for raw in text.splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        p = raw.rstrip("\n").split("\t")
        if len(p) < 12:
            raise AlignmentFormatError(f"Invalid {fmt.upper()} record: {raw[:80]}")
        nmatch, blen = int(p[9]), int(p[10])
        mapq = int(p[11])
        out.append(
            MappingRecord(
                p[0],
                int(p[1]),
                int(p[2]),
                int(p[3]),
                p[4],
                p[5],
                int(p[6]),
                int(p[7]),
                int(p[8]),
                nmatch / blen if blen else None,
                blen,
                None if mapq == 255 else mapq,
            )
        )
    if not out:
        raise AlignmentFormatError(f"No {fmt.upper()} records found.")
    return MappingAlignment(out, fmt, source)


def parse_psl(text: str, source: str = "") -> MappingAlignment:
    out = []
    for raw in text.splitlines():
        p = raw.split("\t")
        if len(p) < 21 or not p[0].strip().isdigit():
            continue  # header lines
        match, mis, rep = int(p[0]), int(p[1]), int(p[2])
        aligned = match + mis + rep
        out.append(
            MappingRecord(
                p[9],
                int(p[10]),
                int(p[11]),
                int(p[12]),
                p[8][0],
                p[13],
                int(p[14]),
                int(p[15]),
                int(p[16]),
                (match + rep) / aligned if aligned else None,
                aligned + int(p[5]) + int(p[7]),
                gap_openings=int(p[4]) + int(p[6]),
            )
        )
    if not out:
        raise AlignmentFormatError("No PSL records found.")
    return MappingAlignment(out, "psl", source)


def parse_blast_tab(text: str, source: str = "") -> MappingAlignment:
    """BLAST/DIAMOND tabular output with the default 12 columns (``-outfmt 6``)."""

    out = []
    for raw in text.splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        p = raw.split("\t")
        if len(p) < 12:
            raise AlignmentFormatError(f"Invalid BLAST tabular record: {raw[:80]}")
        qs, qe, ss, se = int(p[6]), int(p[7]), int(p[8]), int(p[9])
        strand = "-" if (qe - qs) * (se - ss) < 0 else "+"
        out.append(
            MappingRecord(
                p[0],
                None,
                min(qs, qe) - 1,
                max(qs, qe),
                strand,
                p[1],
                None,
                min(ss, se) - 1,
                max(ss, se),
                float(p[2]) / 100,
                int(p[3]),
                evalue=float(p[10]),
                bitscore=float(p[11]),
                gap_openings=int(p[5]),
            )
        )
    if not out:
        raise AlignmentFormatError("No BLAST tabular records found.")
    return MappingAlignment(out, "blast_tab", source)
