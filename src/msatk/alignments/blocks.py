"""Block-structured (whole-genome / comparative) alignments: MAF, XMFA, AXT.

A block alignment is a list of local alignment blocks, each with genomic
coordinates for every row. Blocks are analysed individually and summarised per
genome; they are never concatenated into one matrix, because that would discard
coordinates, strand and block membership.

Coordinates are stored 0-based, half-open, on the forward strand of each source
sequence. Genome names follow the MAF ``genome.chromosome`` convention (text
before the first ``.``); XMFA rows are named by their sequence index or file.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from statistics import median
from typing import Any

from msatk.exceptions import AlignmentFormatError
from msatk.models import Alignment, SequenceRecord


@dataclass
class BlockRow:
    src: str
    start: int  # 0-based forward-strand start
    end: int  # half-open
    strand: str
    src_size: int | None
    text: str

    @property
    def genome(self) -> str:
        return self.src.split(".", 1)[0]


@dataclass
class AlignmentBlock:
    rows: list[BlockRow]
    score: float | None = None

    @property
    def width(self) -> int:
        return max((len(r.text) for r in self.rows), default=0)

    def to_alignment(self) -> Alignment:
        seen: dict[str, int] = defaultdict(int)
        records = []
        for r in self.rows:
            seen[r.src] += 1
            rid = r.src if seen[r.src] == 1 else f"{r.src}#{seen[r.src]}"
            records.append(
                SequenceRecord(rid, r.text.upper(), f"{r.src}:{r.start}-{r.end}{r.strand}")
            )
        return Alignment(tuple(records), fmt="block")


@dataclass
class BlockAlignment:
    kind = "block"
    blocks: list[AlignmentBlock]
    fmt: str
    source: str = ""
    meta: dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------ helpers
    @property
    def genomes(self) -> list[str]:
        order: dict[str, None] = {}
        for b in self.blocks:
            for r in b.rows:
                order.setdefault(r.genome, None)
        return list(order)

    @property
    def reference(self) -> str | None:
        return self.blocks[0].rows[0].genome if self.blocks and self.blocks[0].rows else None

    def block_rows(self, molecule_type: str = "dna") -> list[dict[str, Any]]:
        from msatk.stats.alignment import compute_statistics

        out = []
        for i, b in enumerate(self.blocks, start=1):
            ref = b.rows[0]
            s = compute_statistics(b.to_alignment(), molecule_type, pairwise="none").summary
            out.append(
                {
                    "block": i,
                    "columns": b.width,
                    "rows": len(b.rows),
                    "genomes": len({r.genome for r in b.rows}),
                    "reference_src": ref.src,
                    "reference_start": ref.start,
                    "reference_end": ref.end,
                    "reference_strand": ref.strand,
                    "gap_fraction": s["gap_fraction"],
                    "variable_sites": s["variable_sites"],
                    "mean_consensus_identity": s["mean_consensus_identity"],
                    "score": b.score,
                }
            )
        return out

    def genome_rows(self) -> list[dict[str, Any]]:
        intervals: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
        sizes: dict[str, dict[str, int]] = defaultdict(dict)
        minus: dict[str, int] = defaultdict(int)
        blocks_present: dict[str, int] = defaultdict(int)
        copies: dict[str, int] = defaultdict(int)
        for b in self.blocks:
            per_block: dict[str, int] = defaultdict(int)
            for r in b.rows:
                intervals[r.genome].append((r.src, r.start, r.end))
                if r.src_size:
                    sizes[r.genome][r.src] = r.src_size
                minus[r.genome] += r.strand == "-"
                per_block[r.genome] += 1
            for g, k in per_block.items():
                blocks_present[g] += 1
                copies[g] += k > 1
        ref = self.reference
        out = []
        for g in self.genomes:
            covered = _union_length(intervals[g])
            total_size = sum(sizes[g].values()) if sizes[g] else None
            aligned = sum(e - s for _, s, e in intervals[g])
            out.append(
                {
                    "genome": g,
                    "is_reference": g == ref,
                    "blocks_present": blocks_present[g],
                    "block_fraction": blocks_present[g] / len(self.blocks) if self.blocks else None,
                    "aligned_bases": aligned,
                    "covered_bases": covered,
                    "duplicated_bases": aligned - covered,
                    "genome_size": total_size,
                    "aligned_fraction": covered / total_size if total_size else None,
                    "minus_strand_rows": minus[g],
                    "blocks_with_multiple_copies": copies[g],
                    "collinearity_breakpoints": self._breakpoints(g) if g != ref else 0,
                }
            )
        return out

    def _breakpoints(self, genome: str) -> int:
        """Adjacent blocks (in reference order) whose order, chromosome or strand in
        ``genome`` is inconsistent with the reference: a simple rearrangement count."""

        ref = self.reference
        pairs = []
        for b in self.blocks:
            r0 = next((r for r in b.rows if r.genome == ref), None)
            g = [r for r in b.rows if r.genome == genome]
            if r0 is not None and len(g) == 1:
                pairs.append((r0.src, r0.start, g[0].src, g[0].start, g[0].strand))
        pairs.sort()
        breaks = 0
        for prev, nxt in zip(pairs, pairs[1:]):
            if prev[0] != nxt[0]:
                continue
            same_dir = prev[4] == nxt[4]
            ordered = (nxt[3] > prev[3]) if prev[4] == "+" else (nxt[3] < prev[3])
            if prev[2] != nxt[2] or not same_dir or not ordered:
                breaks += 1
        return breaks

    def summary(self, molecule_type: str = "dna") -> dict[str, Any]:
        widths = sorted(b.width for b in self.blocks)
        genomes = self.genome_rows()
        ref = next((g for g in genomes if g["is_reference"]), None)
        return {
            "alignment_kind": "block",
            "detected_format": self.fmt,
            "blocks": len(self.blocks),
            "genomes": len(genomes),
            "reference_genome": self.reference,
            "total_block_columns": sum(widths),
            "min_block_columns": widths[0] if widths else None,
            "median_block_columns": median(widths) if widths else None,
            "max_block_columns": widths[-1] if widths else None,
            "n50_block_columns": _n50(widths),
            "mean_genomes_per_block": sum(len({r.genome for r in b.rows}) for b in self.blocks)
            / len(self.blocks)
            if self.blocks
            else None,
            "blocks_with_all_genomes": sum(
                1 for b in self.blocks if len({r.genome for r in b.rows}) == len(genomes)
            ),
            "reference_aligned_fraction": ref["aligned_fraction"] if ref else None,
            "minus_strand_rows": sum(g["minus_strand_rows"] for g in genomes),
            "collinearity_breakpoints": sum(g["collinearity_breakpoints"] for g in genomes),
            "duplicated_bases": sum(g["duplicated_bases"] for g in genomes),
        }

    def tables(self, molecule_type: str = "dna") -> dict[str, list[dict[str, Any]]]:
        return {"blocks": self.block_rows(molecule_type), "genomes": self.genome_rows()}


def _union_length(intervals: list[tuple[str, int, int]]) -> int:
    total = 0
    by_src: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for src, s, e in intervals:
        by_src[src].append((s, e))
    for spans in by_src.values():
        spans.sort()
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s > cur_e:
                total += cur_e - cur_s
                cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        total += cur_e - cur_s
    return total


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


# --------------------------------------------------------------------------- parsers
def parse_maf(text: str, source: str = "") -> BlockAlignment:
    blocks: list[AlignmentBlock] = []
    current: AlignmentBlock | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("a"):
            m = re.search(r"score=([-\d.eE+]+)", line)
            current = AlignmentBlock([], float(m.group(1)) if m else None)
            blocks.append(current)
        elif line.startswith("s ") and current is not None:
            parts = line.split()
            if len(parts) < 7:
                raise AlignmentFormatError("Invalid MAF 's' line.")
            src, start, size, strand, src_size, text_ = (
                parts[1],
                int(parts[2]),
                int(parts[3]),
                parts[4],
                int(parts[5]),
                parts[6],
            )
            if strand == "-":  # MAF minus-strand starts are on the reverse complement
                start = src_size - start - size
            current.rows.append(BlockRow(src, start, start + size, strand, src_size, text_))
    blocks = [b for b in blocks if b.rows]
    if not blocks:
        raise AlignmentFormatError("No MAF alignment blocks found.")
    return BlockAlignment(blocks, "maf", source)


def parse_xmfa(text: str, source: str = "") -> BlockAlignment:
    blocks: list[AlignmentBlock] = []
    rows: list[BlockRow] = []
    header: tuple[str, int, int, str] | None = None
    seq: list[str] = []
    names: dict[str, str] = {}

    def flush_row() -> None:
        if header is not None:
            src, s, e, strand = header
            if e > 0:  # 0-0 marks a genome absent from the block
                rows.append(BlockRow(names.get(src, src), s - 1, e, strand, None, "".join(seq)))

    for raw in text.splitlines():
        line = raw.strip()
        m = re.match(r"#Sequence(\d+)File\s+(\S+)", line)
        if m:
            names[m.group(1)] = m.group(2)
            continue
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            flush_row()
            m = re.match(r">\s*(\S+?):(\d+)-(\d+)\s+([+-])", line)
            if not m:
                raise AlignmentFormatError(f"Invalid XMFA header: {line}")
            header = (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
            seq = []
        elif line.startswith("="):
            flush_row()
            header, seq = None, []
            if rows:
                blocks.append(AlignmentBlock(rows))
            rows = []
        else:
            seq.append(line)
    flush_row()
    if rows:
        blocks.append(AlignmentBlock(rows))
    if not blocks:
        raise AlignmentFormatError("No XMFA alignment blocks found.")
    return BlockAlignment(blocks, "xmfa", source)


def parse_axt(text: str, source: str = "") -> BlockAlignment:
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    blocks = []
    for i in range(0, len(lines) - 2, 3):
        h = lines[i].split()
        if len(h) < 9:
            raise AlignmentFormatError(f"Invalid AXT header: {lines[i]}")
        t_row = BlockRow(h[1], int(h[2]) - 1, int(h[3]), "+", None, lines[i + 1])
        q_row = BlockRow(h[4], int(h[5]) - 1, int(h[6]), h[7], None, lines[i + 2])
        blocks.append(AlignmentBlock([t_row, q_row], float(h[8])))
    if not blocks:
        raise AlignmentFormatError("No AXT alignment blocks found.")
    return BlockAlignment(blocks, "axt", source)
