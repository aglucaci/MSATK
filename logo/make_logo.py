"""Draw the MSATK logo (v4) as SVG and PNG, in light and dark variants.

The mark is a real five-sequence alignment beside its guide tree. Each cell is
coloured by consensus support: teal where every sequence agrees, light teal
where most do, grey for substitutions and weakly supported columns,
so the picture shows what MSATK measures. Text is converted to paths, so the
SVG files need no fonts installed.

Usage::

    python logo/make_logo.py          # writes logo/msatk_*.{svg,png}
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.font_manager import FontProperties, findfont  # noqa: E402
from matplotlib.patches import FancyBboxPatch, PathPatch  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402

HERE = Path(__file__).resolve().parent
plt.rcParams["svg.fonttype"] = "path"
plt.rcParams["svg.hashsalt"] = "msatk-logo"  # deterministic SVG ids

ALIGNMENT = [
    "ATGCTA-",
    "ATG-TAG",
    "ACGCTG-",
    "AT-CTAG",
    "ATGCTA-",
]
# ((s1, s2), (s3, s4)), s5 - as in the original v1 logo.
PARTIAL = 0.6  # consensus support for the lighter highlight (full support = 1.0)

THEMES = {
    "light": {
        "ink": "#0B2A5B",
        "soft": "#5B6B85",
        "hit": "#14A39A",
        "part": "#9ED8D1",
        "cell": "#DCE3EE",
    },
    "dark": {
        "ink": "#E8EEF7",
        "soft": "#9FB0C8",
        "hit": "#2DD4BF",
        "part": "#1F6F6A",
        "cell": "#2A3B57",
    },
}

CELL, GAP = 36.0, 9.0
PITCH = CELL + GAP


def _font(family: str, weight: str) -> FontProperties:
    prop = FontProperties(family=family, weight=weight)
    findfont(prop, fallback_to_default=False)
    return prop


def _consensus(columns: list[str]) -> list[tuple[str, float]]:
    out = []
    for col in columns:
        residues = [c for c in col if c != "-"]
        top, n = Counter(residues).most_common(1)[0]
        out.append((top, n / len(col)))
    return out


def _draw_grid(ax, x0: float, y0: float, theme: dict[str, str]) -> None:
    columns = ["".join(row[j] for row in ALIGNMENT) for j in range(len(ALIGNMENT[0]))]
    cons = _consensus(columns)
    for i, row in enumerate(ALIGNMENT):
        y = y0 - i * PITCH
        for j, ch in enumerate(row):
            x = x0 + j * PITCH
            if ch == "-":
                ax.add_patch(
                    FancyBboxPatch(
                        (x + CELL * 0.28, y + CELL / 2 - 2.5),
                        CELL * 0.44,
                        5,
                        boxstyle="round,pad=0,rounding_size=2.5",
                        fc=theme["soft"],
                        ec="none",
                    )
                )
                continue
            top, support = cons[j]
            if ch != top or support < PARTIAL:
                colour = theme["cell"]
            else:
                colour = theme["hit"] if support == 1.0 else theme["part"]
            ax.add_patch(
                FancyBboxPatch(
                    (x, y),
                    CELL,
                    CELL,
                    boxstyle="round,pad=0,rounding_size=8",
                    fc=colour,
                    ec="none",
                )
            )


def _draw_tree(ax, x_right: float, y0: float, theme: dict[str, str]) -> float:
    """Draw ((1,2),(3,4)),5 with leaves ending at x_right. Returns the root x."""
    lw = 5.0
    yc = [y0 - i * PITCH + CELL / 2 for i in range(len(ALIGNMENT))]
    step = 30.0
    x1, x2, x3 = x_right - step, x_right - 2 * step, x_right - 3 * step
    kw = {"color": theme["ink"], "lw": lw, "solid_capstyle": "round", "solid_joinstyle": "round"}
    a = (yc[0] + yc[1]) / 2
    b = (yc[2] + yc[3]) / 2
    c = (a + b) / 2
    root = (c + yc[4]) / 2
    for y in yc[:4]:
        ax.plot([x1, x_right], [y, y], **kw)
    ax.plot([x1, x1], [yc[0], yc[1]], **kw)
    ax.plot([x1, x1], [yc[2], yc[3]], **kw)
    ax.plot([x2, x1], [a, a], **kw)
    ax.plot([x2, x1], [b, b], **kw)
    ax.plot([x2, x2], [a, b], **kw)
    ax.plot([x3, x2], [c, c], **kw)
    ax.plot([x3, x_right], [yc[4], yc[4]], **kw)
    ax.plot([x3, x3], [c, yc[4]], **kw)
    ax.plot([x3 - 14, x3], [root, root], **kw)
    return x3 - 14


def _text(ax, s: str, x: float, y: float, size: float, prop: FontProperties, color: str) -> float:
    path = TextPath((x, y), s, size=size, prop=prop)
    ax.add_patch(PathPatch(path, fc=color, ec="none", lw=0))
    ext = path.get_extents()
    return ext.x1


def _save(fig, ax, stem: str, box: tuple[float, float, float, float], px_width: int) -> None:
    x0, y0, x1, y1 = box
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    w, h = x1 - x0, y1 - y0
    fig.set_size_inches(w / 100, h / 100)
    fig.savefig(HERE / f"{stem}.svg", transparent=True)
    fig.savefig(HERE / f"{stem}.png", transparent=True, dpi=100 * px_width / w)
    plt.close(fig)


def _new_axes():
    fig = plt.figure()
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    ax.set_aspect("equal")
    return fig, ax


def full_logo(variant: str) -> None:
    theme = THEMES[variant]
    fig, ax = _new_axes()
    rows = len(ALIGNMENT)
    top = (rows - 1) * PITCH
    grid_x = 0.0
    _draw_grid(ax, grid_x, top, theme)
    root_x = _draw_tree(ax, grid_x - 14, top, theme)
    grid_right = grid_x + len(ALIGNMENT[0]) * PITCH - GAP
    block_h = rows * PITCH - GAP
    # Divider
    div_x = grid_right + 44
    ax.plot([div_x, div_x], [0, block_h], color=theme["ink"], lw=3, solid_capstyle="round")
    # Wordmark: MSA bold, TK book weight, cap height matched to the block.
    bold = _font("Avenir Next", "bold")
    book = _font("Avenir", "normal")
    size = 190.0
    cap = TextPath((0, 0), "M", size=size, prop=bold).get_extents().height
    baseline = block_h / 2 - cap / 2
    x = _text(ax, "MSA", div_x + 44, baseline, size, bold, theme["ink"])
    x = _text(ax, "TK", x + 6, baseline, size, book, theme["ink"])
    pad = 40
    _save(
        fig,
        ax,
        f"msatk_logo{'' if variant == 'light' else '_dark'}",
        (root_x - pad, -pad, x + pad, block_h + pad),
        px_width=1600,
    )


def mark(variant: str) -> None:
    theme = THEMES[variant]
    fig, ax = _new_axes()
    rows, cols = len(ALIGNMENT), len(ALIGNMENT[0])
    top = (rows - 1) * PITCH
    _draw_grid(ax, 0, top, theme)
    w = cols * PITCH - GAP
    h = rows * PITCH - GAP
    side = max(w, h) + 2 * 24
    cx, cy = w / 2, h / 2
    _save(
        fig,
        ax,
        f"msatk_mark{'' if variant == 'light' else '_dark'}",
        (cx - side / 2, cy - side / 2, cx + side / 2, cy + side / 2),
        px_width=512,
    )


def main() -> None:
    for variant in THEMES:
        full_logo(variant)
        mark(variant)
    for f in sorted(HERE.glob("msatk_*.*")):
        print(f.relative_to(HERE.parent))


if __name__ == "__main__":
    main()
