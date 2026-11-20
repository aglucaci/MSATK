"""Publication-quality residue-level alignment figures.

``plot_alignment`` renders a region of the alignment as a coloured residue grid
with optional letters, sequence labels, a consensus row and conservation /
gap tracks. Output format follows the file suffix (PNG, SVG, PDF, TIFF).
For large alignments render a region (``start``/``end``, ``max_rows``) rather
than compressing the whole matrix into an unreadable image.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from msatk.core.encoded import encode_alignment
from msatk.encoding.images import rgb_image
from msatk.encoding.transforms import subset_alignment
from msatk.exceptions import MissingDependencyError
from msatk.models import Alignment
from msatk.stats.site import site_arrays

SUPPORTED_FORMATS = ("png", "svg", "pdf", "tif", "tiff", "jpg", "jpeg")


def _pyplot() -> Any:
    try:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise MissingDependencyError(
            "Figures require matplotlib: pip install 'msatk[plots]'"
        ) from exc
    return plt


def plot_alignment(
    alignment: Alignment,
    output: str | Path,
    molecule_type: str,
    start: int | None = None,
    end: int | None = None,
    rows: list[str] | None = None,
    max_rows: int | None = 200,
    order: str = "input",
    theme: str | dict[str, str] | None = "auto",
    show_letters: bool | None = None,
    show_consensus: bool = True,
    tracks: tuple[str, ...] = ("conservation", "gaps"),
    highlight_differences: bool = False,
    font_family: str = "DejaVu Sans Mono",
    font_size: float | None = None,
    cell_width: float = 0.14,
    cell_height: float = 0.16,
    dpi: int = 300,
    title: str | None = None,
) -> Path:
    """Render a residue-level alignment figure; returns the written path.

    ``order``: ``input`` | ``consensus_identity`` (most similar to consensus first) |
    ``gaps`` (fewest gaps first). ``highlight_differences`` greys out residues that
    match the consensus so substitutions stand out.
    """

    plt = _pyplot()
    target = Path(output)
    suffix = target.suffix.lower().lstrip(".")
    if suffix not in SUPPORTED_FORMATS:
        raise ValueError(f"Unsupported figure format {suffix!r}; use one of {SUPPORTED_FORMATS}.")

    full_enc = encode_alignment(alignment, molecule_type)
    consensus_full = full_enc.consensus_string()
    order_ids = None
    if order != "input":
        from msatk.stats.sequence import per_sequence_stats_encoded

        rows_stats = per_sequence_stats_encoded(full_enc)
        if order == "consensus_identity":
            key = lambda r: -(r["consensus_identity"] or 0.0)  # noqa: E731
        elif order == "gaps":
            key = lambda r: r["gap_fraction"]  # noqa: E731
        else:
            raise ValueError("order must be input, consensus_identity or gaps")
        order_ids = [r["sequence_id"] for r in sorted(rows_stats, key=key)]
    sub, coords = subset_alignment(alignment, start=start, end=end, ids=rows, order=order_ids)
    if max_rows is not None and sub.n_sequences > max_rows:
        sub, extra = subset_alignment(sub, ids=sub.ids[:max_rows])
        coords.operations.append({"op": "truncate_rows", "max_rows": max_rows})
    enc = encode_alignment(sub, molecule_type)
    n, width = enc.n, enc.width
    img = rgb_image(sub, molecule_type, theme=theme).astype(float) / 255.0
    consensus = "".join(consensus_full[s - 1] for s in coords.original_sites)
    if highlight_differences:
        same = np.array(
            [[ch == consensus[j] for j, ch in enumerate(row)] for row in enc.row_strings()]
        )
        img[same] = img[same] * 0.25 + 0.75
    if show_consensus:
        cons_alignment = Alignment(
            (type(sub.records[0])("consensus", consensus),), source="", fmt="consensus"
        )
        cons_img = rgb_image(cons_alignment, molecule_type, theme=theme).astype(float) / 255.0
        img = np.concatenate([cons_img, np.ones((1, width, 3)), img], axis=0)
    labels = (["Consensus", ""] if show_consensus else []) + sub.ids
    letters = show_letters if show_letters is not None else (width <= 150 and n <= 120)

    track_data = []
    arrays = site_arrays(enc)
    for name in tracks:
        if name == "conservation":
            track_data.append(
                ("Conservation", np.nan_to_num(arrays["dominant_frequency"], nan=0.0))
            )
        elif name == "gaps":
            track_data.append(("Gap fraction", arrays["gaps"] / max(n, 1)))
        elif name == "entropy":
            track_data.append(
                ("Norm. entropy", np.nan_to_num(arrays["normalized_entropy"], nan=0.0))
            )
        elif name == "occupancy":
            track_data.append(("Occupancy", arrays["observed"] / max(n, 1)))
        else:
            raise ValueError(f"Unknown track {name!r}")

    label_w = min(2.6, 0.075 * max(len(lab) for lab in labels) + 0.3)
    fig_w = label_w + max(3.0, width * cell_width) + 0.4
    grid_h = img.shape[0] * cell_height
    track_h = 0.55 * len(track_data)
    fig_h = grid_h + track_h + 0.9
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    heights = [grid_h] + [0.55] * len(track_data)
    gs = fig.add_gridspec(
        len(heights), 1, height_ratios=heights, hspace=0.08,
        left=label_w / fig_w, right=1 - 0.2 / fig_w, top=1 - 0.55 / fig_h, bottom=0.35 / fig_h,
    )  # fmt: skip
    ax = fig.add_subplot(gs[0])
    ax.imshow(
        img,
        aspect="auto",
        interpolation="nearest",
        extent=(0.5, width + 0.5, img.shape[0] - 0.5, -0.5),
    )
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(
        labels,
        fontsize=font_size or max(4.0, min(8.0, 160 / max(len(labels), 1))),
        family=font_family,
    )
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    tick_step = _tick_step(width)
    positions = [
        i for i in range(1, width + 1) if coords.original_sites[i - 1] % tick_step == 0 or i == 1
    ]
    ax.set_xticks(positions)
    ax.set_xticklabels([str(coords.original_sites[i - 1]) for i in positions], fontsize=6)
    ax.xaxis.set_ticks_position("top")
    if letters:
        rows_text = ([consensus, ""] if show_consensus else []) + enc.row_strings()
        fs = font_size or max(3.0, min(7.0, cell_width * 50))
        for r, text in enumerate(rows_text):
            for c, ch in enumerate(text):
                if ch in "-.":
                    continue
                ax.text(
                    c + 1,
                    r,
                    ch,
                    ha="center",
                    va="center",
                    fontsize=fs,
                    family=font_family,
                    color="#111",
                )
    if title:
        ax.set_title(title, fontsize=9, pad=14)
    for k, (name, values) in enumerate(track_data, start=1):
        tax = fig.add_subplot(gs[k], sharex=ax)
        tax.bar(
            np.arange(1, width + 1),
            values,
            width=1.0,
            color="#0f766e" if "Gap" not in name else "#9a3412",
        )
        tax.set_ylim(0, 1)
        tax.set_yticks([0, 1])
        tax.tick_params(labelsize=5, length=2)
        tax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=6)
        tax.tick_params(axis="x", labelbottom=k == len(track_data))
        for side in ("top", "right"):
            tax.spines[side].set_visible(False)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=dpi)
    plt.close(fig)
    return target


def _tick_step(width: int) -> int:
    for step in (1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000):
        if width / step <= 25:
            return step
    return 10000
