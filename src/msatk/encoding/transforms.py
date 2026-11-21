"""Coordinate-preserving alignment transforms and fixed-size windowing.

Biological alignments are not photographs: rotation, mirroring, interpolation
and colour jitter destroy their meaning. MSATK therefore offers only
integer-exact operations (row/column selection, reordering, windowing with
explicit padding masks) and records every operation so that each output
coordinate maps back to the original alignment.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from msatk.models import Alignment, SequenceRecord


@dataclass
class CoordinateMap:
    """Maps rows/sites of a transformed alignment back to the original (1-based sites)."""

    original_rows: list[int]
    original_sites: list[int]
    operations: list[dict[str, Any]] = field(default_factory=list)

    def site_rows(self) -> list[dict[str, int]]:
        return [
            {"output_site": i + 1, "original_site": s} for i, s in enumerate(self.original_sites)
        ]

    def as_dict(self) -> dict[str, Any]:
        return {
            "original_rows_0based": self.original_rows,
            "original_sites_1based": self.original_sites,
            "operations": self.operations,
        }


def subset_alignment(
    alignment: Alignment,
    start: int | None = None,
    end: int | None = None,
    ids: list[str] | None = None,
    order: list[str] | None = None,
    max_site_gap_fraction: float | None = None,
) -> tuple[Alignment, CoordinateMap]:
    """Select a 1-based inclusive site range, a subset/ordering of rows and optionally
    drop gap-rich columns, returning the new alignment and its coordinate map."""

    seqs = alignment.padded_sequences()
    width = alignment.length
    s = 1 if start is None else max(1, start)
    e = width if end is None else min(width, end)
    if s > e:
        raise ValueError(f"Empty site range: start={start} end={end} (width {width}).")
    ops: list[dict[str, Any]] = []
    if (s, e) != (1, width):
        ops.append({"op": "crop_sites", "start": s, "end": e})
    rows = list(range(alignment.n_sequences))
    if ids is not None:
        wanted = set(ids)
        rows = [i for i in rows if alignment.ids[i] in wanted]
        ops.append({"op": "select_rows", "n": len(rows)})
    if order is not None:
        rank = {seq_id: k for k, seq_id in enumerate(order)}
        rows.sort(key=lambda i: rank.get(alignment.ids[i], len(rank) + i))
        ops.append({"op": "reorder_rows"})
    sites = list(range(s, e + 1))
    if max_site_gap_fraction is not None:
        kept = []
        for site in sites:
            col = [seqs[i][site - 1] for i in rows]
            if sum(ch in "-." for ch in col) / max(len(col), 1) <= max_site_gap_fraction:
                kept.append(site)
        ops.append(
            {
                "op": "drop_gappy_sites",
                "max_site_gap_fraction": max_site_gap_fraction,
                "removed": len(sites) - len(kept),
            }
        )
        sites = kept
    if not rows or not sites:
        raise ValueError("Transform removed every row or every site.")
    idx = np.array(sites) - 1
    records = tuple(
        SequenceRecord(
            alignment.records[i].id,
            "".join(np.array(list(seqs[i]))[idx]),
            alignment.records[i].description,
        )
        for i in rows
    )
    return Alignment(records, source=alignment.source, fmt=alignment.fmt), CoordinateMap(
        rows, sites, ops
    )


@dataclass
class Windows:
    """Fixed-size tiles cut from an ``(N, L, ...)`` tensor."""

    data: np.ndarray  # (T, rows, width, ...)
    mask: np.ndarray  # (T, rows, width) True where the cell is real (not padding)
    table: list[dict[str, int]]  # per-tile provenance


def sliding_windows(
    tensor: np.ndarray,
    width: int,
    stride: int | None = None,
    rows: int | None = None,
    row_stride: int | None = None,
    pad_value: float = 0,
    site_offset: int = 0,
) -> Windows:
    """Tile an ``(N, L, ...)`` tensor into ``(rows x width)`` windows without resampling.

    The final window along each axis is padded with ``pad_value`` and the padding
    is marked ``False`` in ``mask``. ``table`` records the original 1-based,
    inclusive site range and 0-based row range of each tile.
    """

    n, length = tensor.shape[:2]
    stride = stride or width
    rows = rows or n
    row_stride = row_stride or rows
    if width < 1 or stride < 1 or rows < 1 or row_stride < 1:
        raise ValueError("Window sizes and strides must be positive integers.")
    col_starts = list(range(0, max(length - width, 0) + 1, stride))
    if col_starts[-1] + width < length:
        col_starts.append(col_starts[-1] + stride)
    row_starts = list(range(0, max(n - rows, 0) + 1, row_stride))
    if row_starts[-1] + rows < n:
        row_starts.append(row_starts[-1] + row_stride)
    tiles: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    table: list[dict[str, int]] = []
    tail = tensor.shape[2:]
    for r0 in row_starts:
        for c0 in col_starts:
            block = np.full((rows, width, *tail), pad_value, dtype=tensor.dtype)
            m = np.zeros((rows, width), dtype=bool)
            sub = tensor[r0 : r0 + rows, c0 : c0 + width]
            block[: sub.shape[0], : sub.shape[1]] = sub
            m[: sub.shape[0], : sub.shape[1]] = True
            tiles.append(block)
            masks.append(m)
            table.append(
                {
                    "tile": len(table),
                    "row_start": r0,
                    "row_end": r0 + sub.shape[0] - 1,
                    "site_start": site_offset + c0 + 1,
                    "site_end": site_offset + c0 + sub.shape[1],
                    "valid_rows": int(sub.shape[0]),
                    "valid_sites": int(sub.shape[1]),
                }
            )
    return Windows(np.stack(tiles), np.stack(masks), table)
