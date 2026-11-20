"""Sequence-graph / pangenome input: GFA 1.x and rGFA.

A variation graph is analysed as a graph — segments, links, paths and walks — and
never flattened into a rectangular matrix, which would lose alternate paths.
Pangenome summaries classify segments by how many paths/walks traverse them:
*core* (all), *shell* (more than one but not all) and *private* (exactly one).
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from statistics import mean
from typing import Any

from msatk.exceptions import AlignmentFormatError


@dataclass
class GraphAlignment:
    kind = "graph"
    segments: dict[str, int]  # name -> length
    links: list[tuple[str, str, str, str]]  # from, from_orient, to, to_orient
    paths: dict[str, list[str]]  # path/walk name -> ordered segment names
    fmt: str = "gfa"
    source: str = ""
    meta: dict[str, Any] = field(default_factory=dict)
    segment_rank: dict[str, int] = field(default_factory=dict)  # rGFA SR tag

    def _degrees(self) -> tuple[dict[tuple[str, str], int], dict[str, int]]:
        ends: dict[tuple[str, str], int] = defaultdict(int)
        deg: dict[str, int] = defaultdict(int)
        for a, ao, b, bo in self.links:
            ends[(a, "R" if ao == "+" else "L")] += 1
            ends[(b, "L" if bo == "+" else "R")] += 1
            deg[a] += 1
            deg[b] += 1
        return ends, deg

    def _components(self) -> int:
        parent = {s: s for s in self.segments}

        def find(x: str) -> str:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a, _, b, _ in self.links:
            if a in parent and b in parent:
                parent[find(a)] = find(b)
        return len({find(s) for s in self.segments})

    def segment_rows(self) -> list[dict[str, Any]]:
        ends, deg = self._degrees()
        support: dict[str, int] = defaultdict(int)
        for segs in self.paths.values():
            for s in set(segs):
                support[s] += 1
        return [
            {
                "segment": s,
                "length": n,
                "degree": deg.get(s, 0),
                "left_links": ends.get((s, "L"), 0),
                "right_links": ends.get((s, "R"), 0),
                "paths": support.get(s, 0),
                "rgfa_rank": self.segment_rank.get(s),
            }
            for s, n in self.segments.items()
        ]

    def path_rows(self) -> list[dict[str, Any]]:
        return [
            {
                "path": name,
                "segments": len(segs),
                "bases": sum(self.segments.get(s, 0) for s in segs),
                "distinct_segments": len(set(segs)),
            }
            for name, segs in self.paths.items()
        ]

    def summary(self) -> dict[str, Any]:
        lengths = sorted(self.segments.values())
        rows = self.segment_rows()
        n_paths = len(self.paths)
        total = sum(lengths)
        core = [r for r in rows if n_paths and r["paths"] == n_paths]
        private = [r for r in rows if r["paths"] == 1]
        shell = [r for r in rows if 1 < r["paths"] < n_paths]
        tips = sum(1 for r in rows if (r["left_links"] == 0) != (r["right_links"] == 0))
        return {
            "alignment_kind": "graph",
            "detected_format": self.fmt,
            "segments": len(lengths),
            "total_segment_bases": total,
            "mean_segment_length": mean(lengths) if lengths else None,
            "n50_segment_length": _n50(lengths),
            "links": len(self.links),
            "connected_components": self._components() if self.segments else 0,
            "isolated_segments": sum(1 for r in rows if r["degree"] == 0),
            "tip_segments": tips,
            "branching_segment_ends": sum(1 for v in self._degrees()[0].values() if v > 1),
            "max_degree": max((r["degree"] for r in rows), default=0),
            "paths": n_paths,
            "mean_path_bases": mean(r["bases"] for r in self.path_rows()) if n_paths else None,
            "segments_on_paths": sum(1 for r in rows if r["paths"]),
            "core_segments": len(core) if n_paths else None,
            "shell_segments": len(shell) if n_paths else None,
            "private_segments": len(private) if n_paths else None,
            "core_bases_fraction": sum(r["length"] for r in core) / total
            if n_paths and total
            else None,
            "reference_rank0_segments": sum(1 for v in self.segment_rank.values() if v == 0)
            if self.segment_rank
            else None,
        }

    def tables(self) -> dict[str, list[dict[str, Any]]]:
        return {"segments": self.segment_rows(), "paths": self.path_rows()}


def _n50(values: list[int]) -> int | None:
    if not values:
        return None
    half, acc = sum(values) / 2, 0
    for v in sorted(values, reverse=True):
        acc += v
        if acc >= half:
            return v
    return None


def parse_gfa(text: str, source: str = "") -> GraphAlignment:
    segments: dict[str, int] = {}
    ranks: dict[str, int] = {}
    links: list[tuple[str, str, str, str]] = []
    paths: dict[str, list[str]] = {}
    version = None
    for raw in text.splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        p = raw.rstrip("\n").split("\t")
        tag = p[0]
        if tag == "H":
            version = next((f[5:] for f in p[1:] if f.startswith("VN:Z:")), None)
        elif tag == "S" and len(p) >= 3:
            ln = next((int(f[5:]) for f in p[3:] if f.startswith("LN:i:")), None)
            segments[p[1]] = ln if ln is not None else (0 if p[2] == "*" else len(p[2]))
            sr = next((int(f[5:]) for f in p[3:] if f.startswith("SR:i:")), None)
            if sr is not None:
                ranks[p[1]] = sr
        elif tag == "L" and len(p) >= 5:
            links.append((p[1], p[2], p[3], p[4]))
        elif tag == "P" and len(p) >= 3:
            paths[p[1]] = [s[:-1] for s in p[2].split(",") if s]
        elif tag == "W" and len(p) >= 7:
            name = f"{p[1]}#{p[2]}#{p[3]}:{p[4]}-{p[5]}"
            paths[name] = re.findall(r"[<>]([^<>]+)", p[6])
    if not segments:
        raise AlignmentFormatError("No GFA segments (S lines) found.")
    if version and version.startswith("2"):
        raise AlignmentFormatError("GFA 2 is not supported yet; convert to GFA 1.")
    return GraphAlignment(segments, links, paths, "gfa", source, {"gfa_version": version}, ranks)
