"""Compose Figure 1 of the Application Note from real MSATK outputs.

Run ``use_cases.py`` first. Writes ``paper/figures/figure1.{png,pdf}``.

Panels
------
A  Architecture: one encoded alignment drives statistics, QC, figures, the
   report and ML encodings, all with provenance.
B  Statistical outputs. Top: per-site information content and gap fraction for
   the Ras family (PF00071) with G-box motifs shaded. Bottom: per-sequence
   identity to consensus for 400 H3N2 HA sequences coloured by host label,
   with QC-flagged divergent sequences outlined.
C  Residue-level alignment figure (sites 1-80, 30 sequences ordered by identity
   to consensus) and a sequence logo of sites 1-20 (P-loop at 6-13).
D  ML representations: lossless label image of the full Ras alignment with the
   64-site x 32-row tiling overlaid, one tile's one-hot tensor, and the manifest
   fields that map every tile back to original coordinates.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

from msatk import MSATK  # noqa: E402

HERE = Path(__file__).parent
RES = HERE / "results"
PANELS = RES / "panels"
OUT = HERE.parent / "figures"

INK = "#1f2933"
MUTED = "#6b7785"
ACCENT = "#0f766e"
GAP = "#9a3412"
HOST_COLOURS = {
    "human": "#9aa5b1",
    "swine": "#d97706",
    "canine": "#7c3aed",
    "environment": "#0891b2",
    "avian": "#16a34a",
}


def label(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.02,
        1.02,
        letter,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="bottom",
        ha="right",
        color=INK,
    )


def _clean(ax: plt.Axes) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=6, length=2)


def panel_a(ax: plt.Axes) -> None:
    ax.set_xlim(-0.1, 20.6)
    ax.set_ylim(-0.1, 5.1)
    ax.axis("off")

    def box(x, y, w, h, text, fc="#f1f5f9", ec="#cbd2d9", size=6.0, weight="normal"):
        ax.add_patch(
            FancyBboxPatch(
                (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.15", fc=fc, ec=ec, lw=0.7
            )
        )
        ax.text(
            x + w / 2,
            y + h / 2,
            text,
            ha="center",
            va="center",
            fontsize=size,
            color=INK,
            fontweight=weight,
            linespacing=1.2,
        )

    def arrow(x0, y0, x1, y1):
        ax.annotate(
            "",
            (x1, y1),
            (x0, y0),
            arrowprops={
                "arrowstyle": "-|>",
                "color": MUTED,
                "lw": 0.8,
                "shrinkA": 0,
                "shrinkB": 0,
                "mutation_scale": 7,
            },
        )

    # Row 1: input -> validation -> encoded alignment
    box(0.0, 3.7, 3.6, 1.3, "Input\nFASTA, PHYLIP, Clustal,\nStockholm, NEXUS")
    box(
        4.4,
        3.7,
        4.6,
        1.3,
        "Validation and detection\nmolecule type, codon frame,\nstops, duplicates, padding",
    )
    box(
        9.8,
        3.7,
        5.2,
        1.3,
        "EncodedAlignment (built once)\nper-cell state: residue, ambiguous,\nmissing, gap, stop, invalid",
        fc="#e6f4f1",
        ec=ACCENT,
    )
    box(
        15.8,
        3.7,
        4.2,
        1.3,
        "Metric registry\n145 definitions, formulae,\nmissing-data rules",
        fc="#fff7ed",
        ec="#f59e0b",
    )
    arrow(3.6, 4.35, 4.4, 4.35)
    arrow(9.0, 4.35, 9.8, 4.35)
    # Row 2: engines
    engines = [
        "Site\nentropy, IC,\nparsimony",
        "Sequence\ngaps, identity,\nrobust z",
        "Pairwise\nexact or sampled\nidentity",
        "Alignment\ncomposition,\ngap blocks",
        "Diversity\nπ, θW, Tajima's D,\nTs/Tv",
        "Codon / protein\nRSCU, frames,\nphysicochemical",
    ]
    w = 3.1
    for i, txt in enumerate(engines):
        box(0.0 + i * (w + 0.28), 1.75, w, 1.35, txt, fc="#ffffff", ec=ACCENT, size=5.8)
    ax.text(
        10.0,
        3.33,
        "Statistics engine (vectorized, shared by every output)",
        ha="center",
        va="center",
        fontsize=6.2,
        color=ACCENT,
        fontweight="bold",
    )
    arrow(12.4, 3.7, 12.4, 3.48)
    # Row 3: outputs
    outs = [
        "Tables: CSV, TSV, JSON",
        "QC findings with items",
        "Figures: PNG, SVG, PDF, TIFF",
        "Offline HTML report",
        "ML tensors + manifest",
    ]
    w3 = 3.75
    for i, txt in enumerate(outs):
        box(0.0 + i * (w3 + 0.31), 0.45, w3, 0.8, txt, fc="#f8fafc", size=5.9)
    arrow(10.0, 1.75, 10.0, 1.25)
    ax.text(
        10.0,
        0.08,
        "Provenance in every run: input SHA-256, versions, parameters, seed (output schema 2.0)",
        ha="center",
        va="center",
        fontsize=5.8,
        color=MUTED,
        style="italic",
    )
    label(ax, "A")


def panel_b(fig: plt.Figure, spec) -> None:
    sub = spec.subgridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.28)
    left = sub[0].subgridspec(2, 1, height_ratios=[2.2, 0.8], hspace=0.12)
    ras = MSATK(HERE / "data" / "PF00071.fasta", sequence_type="protein", pairwise="none").stats()
    uc = json.loads((RES / "use_cases.json").read_text())
    motifs = uc["ras_pf00071"]["g_box_motifs"]
    site = np.array([s["site_index"] for s in ras.sites])
    ic = np.array([s["information_content"] or 0 for s in ras.sites], dtype=float)
    gap = np.array([s["gap_fraction"] for s in ras.sites], dtype=float)
    occupied = gap <= 0.5

    ax = fig.add_subplot(left[0])
    for name, m in motifs.items():
        a, b = m["sites"]
        ax.axvspan(a - 0.5, b + 0.5, color="#fde68a", alpha=0.8, lw=0)
        ax.text(
            (a + b) / 2, 4.45, name.split(" ")[0], ha="center", va="bottom", fontsize=5.8, color=INK
        )
    ax.bar(site[occupied], ic[occupied], width=1.0, color=ACCENT, lw=0, label="gaps ≤ 50%")
    ax.bar(site[~occupied], ic[~occupied], width=1.0, color="#d5dbe1", lw=0, label="gaps > 50%")
    ax.set_xlim(0.5, site[-1] + 0.5)
    ax.set_ylim(0, 4.9)
    ax.set_ylabel("information (bits)", fontsize=6.5)
    ax.set_title("Ras, PF00071 seed (60 × 229)", fontsize=7, loc="left", color=INK)
    ax.legend(
        fontsize=5.5,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(1.02, 1.19),
        ncol=2,
        handlelength=1,
        columnspacing=0.8,
    )
    ax.tick_params(labelbottom=False)
    _clean(ax)
    label(ax, "B")
    ax2 = fig.add_subplot(left[1], sharex=ax)
    ax2.bar(site, gap, width=1.0, color=GAP, lw=0)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([0, 1])
    ax2.set_ylabel("gaps", fontsize=6.5)
    ax2.set_xlabel("alignment site", fontsize=6.5)
    _clean(ax2)

    ax3 = fig.add_subplot(sub[1])
    with (HERE / "data" / "h3n2_ha_metadata.tsv").open() as h:
        host = {r["accession"]: r["host_label"] for r in csv.DictReader(h, delimiter="\t")}
    rows = list(csv.DictReader((RES / "h3n2" / "tables" / "per_sequence_stats.csv").open()))
    findings = json.loads((RES / "h3n2" / "qc_findings.json").read_text())
    findings = findings if isinstance(findings, list) else findings["findings"]
    divergent = set(next(f for f in findings if f["code"] == "DIVERGENT_SEQUENCES")["items"])

    def group(h: str) -> str:
        return h if h in HOST_COLOURS else "avian"

    rows.sort(key=lambda r: float(r["consensus_identity"]))
    x = np.arange(1, len(rows) + 1)
    y = np.array([float(r["consensus_identity"]) for r in rows])
    colours = np.array([HOST_COLOURS[group(host[r["sequence_id"]])] for r in rows])
    flagged = np.array([r["sequence_id"] in divergent for r in rows])
    order = np.argsort([c != HOST_COLOURS["human"] for c in colours])
    nf = order[~flagged[order]]
    ax3.scatter(x[nf], y[nf], s=5, c=colours[nf], lw=0, zorder=2)
    ax3.scatter(
        x[flagged], y[flagged], s=14, c=colours[flagged], edgecolors=INK, linewidths=0.5, zorder=3
    )
    for name, c in HOST_COLOURS.items():
        ax3.scatter([], [], s=10, c=c, label=name)
    ax3.scatter(
        [], [], s=14, facecolors="none", edgecolors=INK, linewidths=0.5, label="QC: divergent"
    )
    ax3.legend(
        fontsize=5.3,
        frameon=False,
        loc="lower right",
        ncol=1,
        handletextpad=0.1,
        labelspacing=0.25,
        borderaxespad=0.1,
    )
    ax3.set_xscale("log")
    ax3.set_ylabel("identity to consensus", fontsize=6.5)
    ax3.set_xlabel("sequence rank (log)", fontsize=6.5)
    ax3.set_title("Influenza A(H3N2) HA CDS (400 × 1,701 nt)", fontsize=7, loc="left", color=INK)
    _clean(ax3)


def panel_c(fig: plt.Figure, spec) -> None:
    sub = spec.subgridspec(1, 2, width_ratios=[1.7, 1.0], wspace=0.04)
    ax = fig.add_subplot(sub[0])
    ax.imshow(plt.imread(PANELS / "ras_alignment.png"), interpolation="lanczos")
    ax.axis("off")
    ax.set_title(
        "Alignment figure: sites 1–50, 20 sequences, conservation and gap tracks",
        fontsize=7,
        loc="left",
        color=INK,
    )
    label(ax, "C")
    ax2 = fig.add_subplot(sub[1])
    ax2.imshow(plt.imread(PANELS / "ras_logo_g1.png"), interpolation="lanczos")
    ax2.axis("off")
    ax2.set_title("Logo, sites 1–20 (P-loop at 6–13)", fontsize=7, loc="left", color=INK)


def panel_d(fig: plt.Figure, spec) -> None:
    sub = spec.subgridspec(1, 3, width_ratios=[1.45, 0.85, 0.95], wspace=0.32)
    data = np.load(RES / "ras_ml" / "ras_onehot.npz")
    manifest = json.loads((RES / "ras_ml" / "ras_onehot.manifest.json").read_text())
    rgb = plt.imread(PANELS / "ras_rgb.png")
    win = manifest["windowing"]

    ax = fig.add_subplot(sub[0])
    ax.imshow(rgb, interpolation="nearest", aspect="auto")
    scale = rgb.shape[1] / 229
    for t in win["tiles"]:
        if t["row_start"] != 0:
            continue
        x0 = (t["site_start"] - 1) * scale
        ax.plot([x0, x0], [0, 32 * scale], color="white", lw=0.7)
    ax.axhline(32 * scale, color="white", lw=0.9)
    ax.add_patch(Rectangle((0, 0), 64 * scale, 32 * scale, fill=False, ec=INK, lw=1.4))
    ax.text(
        32 * scale,
        16 * scale,
        "tile 0",
        color=INK,
        fontsize=6,
        ha="center",
        va="center",
        fontweight="bold",
        bbox={"fc": "white", "ec": "none", "pad": 1, "alpha": 0.85},
    )
    ax.set_xticks(
        [(s - 1) * scale for s in (1, 50, 100, 150, 200)], ["1", "50", "100", "150", "200"]
    )
    ax.set_yticks([])
    ax.set_xlabel("alignment site", fontsize=6.5)
    ax.set_title(
        f"{win['n_tiles']} windows, {win['width']} sites × {win['rows_per_tile']} rows, stride {win['stride']}",
        fontsize=7,
        loc="left",
        color=INK,
    )
    _clean(ax)
    label(ax, "D")

    ax2 = fig.add_subplot(sub[1])
    ch = manifest["channels"]
    ax2.imshow(
        data["x"][0][0].T, aspect="auto", cmap="Greys", interpolation="nearest", vmin=0, vmax=1
    )
    ax2.set_yticks(range(len(ch)), [c.strip("<>") for c in ch], fontsize=4.3)
    ax2.set_xticks([0, 31, 63], ["1", "32", "64"])
    ax2.set_xlabel("site (tile 0, row 1)", fontsize=6.5)
    ax2.set_title("one-hot (22 ch.)", fontsize=7, loc="left", color=INK)
    _clean(ax2)

    ax3 = fig.add_subplot(sub[2])
    ax3.axis("off")
    t6 = win["tiles"][6]
    lines = [
        '"schema": "msatk.encoding/1.0"',
        f'"shape": {manifest["shape"]}',
        '"axes": [tile, seq, site, ch]',
        f'"lossless": {str(manifest["lossless"]).lower()}',
        f'"ambiguity": "{manifest["encoding_details"]["ambiguity_policy"]}"',
        f'"tiles"[6]: sites {t6["site_start"]}–{t6["site_end"]}',
        f"   valid {t6['valid_sites']}/{win['width']} → mask",
        '"input_sha256": "870ee1ec…"',
    ]
    ax3.set_title("manifest.json (excerpt)", fontsize=7, loc="left", color=INK)
    ax3.text(
        0.0,
        0.95,
        "\n".join(lines),
        fontsize=5.2,
        family="DejaVu Sans Mono",
        color=INK,
        va="top",
        transform=ax3.transAxes,
        linespacing=1.5,
        bbox={"boxstyle": "round,pad=0.45", "fc": "#f8fafc", "ec": "#cbd2d9", "lw": 0.6},
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(7.2, 8.9), dpi=300)
    grid = fig.add_gridspec(
        4,
        1,
        height_ratios=[0.95, 1.15, 1.1, 1.0],
        hspace=0.32,
        left=0.075,
        right=0.985,
        top=0.975,
        bottom=0.045,
    )
    panel_a(fig.add_subplot(grid[0]))
    panel_b(fig, grid[1])
    panel_c(fig, grid[2])
    panel_d(fig, grid[3])
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"figure1.{ext}", dpi=300)
    print(f"wrote {OUT / 'figure1.png'}")


if __name__ == "__main__":
    main()
