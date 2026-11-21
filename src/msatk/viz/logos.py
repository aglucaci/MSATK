"""Sequence logos (information-content stacks) rendered with matplotlib."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from msatk.core.encoded import encode_alignment
from msatk.encoding.transforms import subset_alignment
from msatk.models import Alignment
from msatk.viz.alignment import _pyplot
from msatk.viz.themes import load_theme


def logo_matrix(
    alignment: Alignment, molecule_type: str, small_sample_correction: bool = True
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Return (alphabet, letter heights ``(L, K)`` in bits, information ``(L,)``).

    Information content is ``log2 K - (H + e_n)``, where ``e_n = (K-1)/(2 ln2 n)``
    is the Schneider et al. (1986) small-sample correction (optional). Gaps and
    missing data are excluded; columns with no residues get zero height.
    """

    enc = encode_alignment(alignment, molecule_type)
    counts = enc.site_counts.astype(float)
    totals = counts.sum(axis=1, keepdims=True)
    k = enc.k
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(totals > 0, counts / totals, 0.0)
        h = -np.where(p > 0, p * np.log2(np.where(p > 0, p, 1)), 0).sum(axis=1)
        correction = (
            (k - 1) / (2 * math.log(2) * np.maximum(totals[:, 0], 1))
            if small_sample_correction
            else 0.0
        )
    info = np.clip(math.log2(k) - (h + correction), 0, None)
    info[totals[:, 0] == 0] = 0
    return list(enc.policy.canonical), p * info[:, None], info


def plot_logo(
    alignment: Alignment,
    output: str | Path,
    molecule_type: str,
    start: int | None = None,
    end: int | None = None,
    theme: str | dict[str, str] | None = "auto",
    width_per_site: float = 0.22,
    height: float = 1.8,
    dpi: int = 300,
    small_sample_correction: bool = True,
    title: str | None = None,
) -> Path:
    """Write a sequence logo of a region (format from suffix)."""

    plt = _pyplot()
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import PathPatch
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D

    sub, coords = subset_alignment(alignment, start=start, end=end)
    alphabet, heights, _ = logo_matrix(sub, molecule_type, small_sample_correction)
    colors = load_theme(theme, molecule_type)
    width = heights.shape[0]
    fig, ax = plt.subplots(figsize=(max(3.0, width * width_per_site + 0.8), height), dpi=dpi)
    font = FontProperties(family="DejaVu Sans", weight="bold")
    for col in range(width):
        y = 0.0
        for idx in np.argsort(heights[col]):
            hgt = heights[col, idx]
            if hgt <= 1e-6:
                continue
            letter = alphabet[idx]
            tp = TextPath((0, 0), letter, size=1, prop=font)
            bb = tp.get_extents()
            transform = (
                Affine2D()
                .translate(-bb.x0, -bb.y0)
                .scale(0.9 / bb.width, hgt / bb.height)
                .translate(col + 1 - 0.45, y)
            )
            ax.add_patch(
                PathPatch(transform.transform_path(tp), color=colors.get(letter, "#444"), lw=0)
            )
            y += hgt
    ax.set_xlim(0.5, width + 0.5)
    ax.set_ylim(0, math.log2(len(alphabet)))
    ax.set_ylabel("bits", fontsize=7)
    step = max(1, width // 20)
    ticks = list(range(1, width + 1, step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(coords.original_sites[t - 1]) for t in ticks], fontsize=6)
    ax.tick_params(axis="y", labelsize=6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if title:
        ax.set_title(title, fontsize=8)
    fig.tight_layout()
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=dpi)
    plt.close(fig)
    return target
