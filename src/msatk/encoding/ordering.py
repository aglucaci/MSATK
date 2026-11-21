"""Row orderings for exported alignments.

The order of sequences changes the spatial patterns a convolutional or vision
model sees, so it is a modelling choice. MSATK makes it explicit and records the
method, its parameters and the resulting order in every export manifest.

Methods
-------
``input``               original file order.
``reference``           decreasing identity to a reference sequence (default: the first).
``consensus_identity``  decreasing identity to the column consensus.
``gaps``                increasing gap fraction.
``cluster``             leaf order of average-linkage (UPGMA) clustering on p-distance.
a list of IDs, or ``file:PATH`` with one ID per line
                        user-defined; IDs not listed keep their input order at the end.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from msatk.core.encoded import encode_alignment
from msatk.models import Alignment

ORDER_METHODS = ("input", "reference", "consensus_identity", "gaps", "cluster")


def row_order(
    alignment: Alignment,
    molecule_type: str,
    method: str | Sequence[str] = "input",
    reference: str | None = None,
) -> tuple[list[str], dict[str, Any]]:
    """Return ``(ordered_ids, details)`` for ``method`` (see module docstring)."""

    ids = alignment.ids
    if not isinstance(method, str):
        return _user_order(ids, list(method), source="list")
    if method.startswith("file:"):
        path = Path(method[5:])
        listed = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines()]
        return _user_order(ids, [x for x in listed if x and not x.startswith("#")], str(path))
    if method == "input":
        return list(ids), {"method": "input"}
    enc = encode_alignment(alignment, molecule_type)
    if method in ("consensus_identity", "gaps"):
        from msatk.stats.sequence import per_sequence_stats_encoded

        rows = per_sequence_stats_encoded(enc)
        if method == "consensus_identity":
            rows.sort(key=lambda r: -(r["consensus_identity"] or 0.0))
        else:
            rows.sort(key=lambda r: r["gap_fraction"])
        return [str(r["sequence_id"]) for r in rows], {"method": method, "stable_sort": True}
    if method not in ("reference", "cluster"):
        raise ValueError(
            f"Unknown row order {method!r}; choose from {ORDER_METHODS}, a list of IDs "
            "or file:PATH."
        )
    from msatk.stats.pairwise import exact_pairwise

    identity = exact_pairwise(enc).identity
    assert identity is not None
    if method == "reference":
        ref = reference if reference is not None else ids[0]
        if ref not in ids:
            raise ValueError(f"Reference sequence {ref!r} is not in the alignment.")
        r = ids.index(ref)
        ident = np.nan_to_num(identity[r], nan=-1.0)
        ident[r] = 2.0  # reference first
        order = sorted(range(len(ids)), key=lambda i: (-ident[i], i))
        return [ids[i] for i in order], {
            "method": "reference",
            "reference": ref,
            "metric": "identity over sites where both sequences have a residue",
            "undefined_identity": "placed last",
        }
    # cluster
    dist = 1.0 - identity
    undefined = int(np.isnan(dist[np.triu_indices(len(ids), 1)]).sum())
    dist = np.where(np.isnan(dist), 1.0, dist)
    np.fill_diagonal(dist, 0.0)
    leaves = _upgma_leaf_order(dist)
    return [ids[i] for i in leaves], {
        "method": "cluster",
        "linkage": "average (UPGMA)",
        "distance": "p-distance over sites where both sequences have a residue",
        "undefined_distances_set_to": 1.0,
        "pairs_with_undefined_distance": undefined,
    }


def _user_order(ids: list[str], listed: list[str], source: str) -> tuple[list[str], dict[str, Any]]:
    known = set(ids)
    unknown = [x for x in listed if x not in known]
    if unknown:
        raise ValueError(
            f"Row order lists {len(unknown)} ID(s) not in the alignment, e.g. {unknown[:3]}."
        )
    head = list(dict.fromkeys(listed))  # de-duplicate, keep first occurrence
    seen = set(head)
    tail = [x for x in ids if x not in seen]
    return head + tail, {
        "method": "user",
        "source": source,
        "listed": len(head),
        "appended_in_input_order": len(tail),
    }


def _upgma_leaf_order(dist: np.ndarray) -> list[int]:
    """Leaf order of average-linkage clustering.

    Deterministic: ties merge the lowest index pair, and a merged cluster lists
    the lower-index cluster's leaves first. Pure NumPy (no SciPy), so the order
    is identical in every environment. O(N^3) time; fine for a few thousand rows.
    """

    n = dist.shape[0]
    if n <= 2:
        return list(range(n))
    d = dist.astype(np.float64).copy()
    np.fill_diagonal(d, np.inf)
    size = np.ones(n)
    members: list[list[int]] = [[i] for i in range(n)]
    active = np.ones(n, dtype=bool)
    for _ in range(n - 1):
        masked = np.where(active[:, None] & active[None, :], d, np.inf)
        flat = int(np.argmin(masked))
        i, j = divmod(flat, n)
        i, j = min(i, j), max(i, j)
        merged = (d[i] * size[i] + d[j] * size[j]) / (size[i] + size[j])
        d[i, :] = merged
        d[:, i] = merged
        d[i, i] = np.inf
        active[j] = False
        d[j, :] = np.inf
        d[:, j] = np.inf
        members[i] = members[i] + members[j]
        size[i] += size[j]
    return members[int(np.flatnonzero(active)[0])]
