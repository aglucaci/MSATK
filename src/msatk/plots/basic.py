"""Standard statistical figures written in every requested format (PNG, SVG, PDF, TIFF)."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

BRAND = "#0f766e"
ALERT = "#9a3412"
VALID_PLOT_FORMATS = ("png", "svg", "pdf", "tif", "tiff")


def parse_plot_formats(spec: str | Sequence[str] | None) -> list[str]:
    """Parse ``"png,svg"``-style specifications and reject unsupported formats."""

    if spec is None:
        return ["png"]
    items = spec.split(",") if isinstance(spec, str) else list(spec)
    formats = [item.strip().lower().lstrip(".") for item in items if item.strip()]
    if formats == ["none"]:
        return []
    bad = [f for f in formats if f not in VALID_PLOT_FORMATS]
    if bad:
        raise ValueError(f"Unsupported plot format(s): {', '.join(bad)}. Use {VALID_PLOT_FORMATS}.")
    return list(dict.fromkeys(formats))


def _num(values: list[Any]) -> list[float]:
    return [float("nan") if v is None else float(v) for v in values]


def generate_standard_plots(
    out_dir: str | Path,
    per_sequence: list[dict[str, Any]],
    per_site: list[dict[str, Any]],
    composition: list[dict[str, Any]],
    identities: list[list[float | None]] | None = None,
    embeddings: list[dict[str, Any]] | None = None,
    formats: Sequence[str] = ("png",),
    pair_identities: Sequence[float] | None = None,
) -> list[str]:
    """Write the standard figure set. Returns figure stems written in the first format."""

    target = Path(out_dir)
    target.mkdir(parents=True, exist_ok=True)
    generated: list[str] = []
    formats = list(formats) or ["png"]
    try:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt
    except Exception:
        (target / "PLOTS_NOT_GENERATED.txt").write_text(
            "Install MSATK with the plots extra to generate figures: pip install 'msatk[plots]'.\n",
            encoding="utf-8",
        )
        return generated

    def save(stem: str, draw: Callable[[], None]) -> None:
        plt.figure(figsize=(9, 4.8), dpi=150)
        draw()
        plt.tight_layout()
        for fmt in formats:
            plt.savefig(target / f"{stem}.{fmt}", dpi=300 if fmt in {"tif", "tiff"} else 150)
        plt.close()
        generated.append(f"{stem}.{formats[0]}")

    sites = [int(row["site_index"]) for row in per_site]
    save(
        "missingness_by_sequence",
        lambda: _bar(
            plt,
            [str(row["sequence_id"]) for row in per_sequence],
            _num([row["gap_fraction"] for row in per_sequence]),
            "Sequence",
            "Gap fraction",
            "Gap fraction by sequence",
        ),
    )
    save(
        "gap_fraction_by_site",
        lambda: _line(
            plt,
            sites,
            _num([r["gap_fraction"] for r in per_site]),
            "Site",
            "Gap fraction",
            "Gap fraction by site",
            ALERT,
        ),
    )
    save(
        "entropy_by_site",
        lambda: _line(
            plt,
            sites,
            _num([r["entropy"] for r in per_site]),
            "Site",
            "Shannon entropy (bits)",
            "Entropy by site (unobserved sites omitted)",
        ),
    )
    save(
        "conservation_by_site",
        lambda: _line(
            plt,
            sites,
            _num([r["conservation_score"] for r in per_site]),
            "Site",
            "Dominant residue frequency",
            "Conservation by site",
        ),
    )
    save(
        "composition_summary",
        lambda: _bar(
            plt,
            [str(r.get("character", r.get("amino_acid", r.get("codon", "")))) for r in composition],
            [float(r.get("fraction", r.get("frequency", 0.0))) for r in composition],
            "Character",
            "Fraction of cells",
            "Character composition",
        ),
    )
    if identities is not None and len(identities) > 1:
        save("pairwise_identity_heatmap", lambda: _heatmap(plt, identities, "Pairwise identity"))
    raw: list[float | None] = list(pair_identities) if pair_identities is not None else []
    if not raw and identities is not None:
        raw = [
            identities[i][j]
            for i in range(len(identities))
            for j in range(i + 1, len(identities))
            if identities[i][j] is not None
        ]
    values = [float(v) for v in raw if v is not None and v == v]
    if values:
        save(
            "pairwise_identity_distribution",
            lambda: _hist(
                plt, values, "Pairwise identity", "Pairs", "Pairwise identity distribution"
            ),
        )
        save(
            "pairwise_distance_distribution",
            lambda: _hist(
                plt,
                [1 - v for v in values],
                "p-distance",
                "Pairs",
                "Pairwise distance distribution",
            ),
        )
    if embeddings:
        method = str(embeddings[0].get("method", "embedding")).upper()
        save(
            "sequence_embedding",
            lambda: _scatter(
                plt,
                [float(r["embedding_1"]) for r in embeddings],
                [float(r["embedding_2"]) for r in embeddings],
                f"{method} 1",
                f"{method} 2",
                f"Sequence embedding ({method})",
            ),
        )
    return generated


def plot_single_series(
    path: str | Path, x: list[int], y: list[Any], xlabel: str, ylabel: str, title: str
) -> bool:
    """Write a single line plot (format from suffix). Returns False when matplotlib is unavailable."""

    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 4.8), dpi=150)
    _line(plt, x, _num(y), xlabel, ylabel, title)
    plt.tight_layout()
    plt.savefig(target)
    plt.close()
    return True


def plot_single_bar(
    path: str | Path, labels: list[str], values: list[float], xlabel: str, ylabel: str, title: str
) -> bool:
    """Write a single bar plot (format from suffix). Returns False when matplotlib is unavailable."""

    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 4.8), dpi=150)
    _bar(plt, labels, values, xlabel, ylabel, title)
    plt.tight_layout()
    plt.savefig(target)
    plt.close()
    return True


def _style(plt) -> None:
    ax = plt.gca()
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _bar(plt, labels: list[str], values: list[float], xlabel: str, ylabel: str, title: str) -> None:
    plt.bar(range(len(values)), values, color=BRAND)
    if len(labels) <= 40:
        plt.xticks(range(len(labels)), labels, rotation=60, ha="right", fontsize=7)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    _style(plt)


def _line(
    plt, x: list[int], y: list[float], xlabel: str, ylabel: str, title: str, color: str = BRAND
) -> None:
    if len(x) > 400:
        plt.plot(x, y, color=color, linewidth=0.8)
    else:
        plt.plot(x, y, color=color, linewidth=1.4, marker="o", markersize=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    _style(plt)


def _heatmap(plt, matrix: list[list[float | None]], title: str) -> None:
    import numpy as np

    data = np.array([[np.nan if v is None else v for v in row] for row in matrix], dtype=float)
    plt.imshow(data, cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
    plt.colorbar(label="Identity (blank = no comparable sites)")
    plt.title(title)


def _hist(plt, values: list[float], xlabel: str, ylabel: str, title: str) -> None:
    plt.hist(values, bins=min(40, max(5, len(values) // 5)), color=BRAND, edgecolor="white")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    _style(plt)


def _scatter(plt, x: list[float], y: list[float], xlabel: str, ylabel: str, title: str) -> None:
    plt.scatter(x, y, color=BRAND, s=24, alpha=0.8)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    _style(plt)
