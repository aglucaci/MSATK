"""Publication-oriented plots with graceful fallbacks."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path


def generate_standard_plots(
    out_dir: str | Path,
    per_sequence: list[dict[str, object]],
    per_site: list[dict[str, object]],
    composition: list[dict[str, object]],
    identities: list[list[float]],
    embeddings: list[dict[str, object]] | None = None,
) -> list[str]:
    target = Path(out_dir)
    target.mkdir(parents=True, exist_ok=True)
    generated: list[str] = []
    try:
        import matplotlib.pyplot as plt
    except Exception:
        (target / "PLOTS_NOT_GENERATED.txt").write_text(
            "Install MSATK with the plots extra to generate PNG figures: pip install 'msatk[plots]'.\n",
            encoding="utf-8",
        )
        return generated

    def save(name: str, draw: Callable[[], None]) -> None:
        plt.figure(figsize=(9, 4.8), dpi=150)
        draw()
        plt.tight_layout()
        plt.savefig(target / name)
        plt.close()
        generated.append(name)

    save(
        "missingness_by_sequence.png",
        lambda: _bar(
            plt,
            [str(row["sequence_id"]) for row in per_sequence],
            [float(row["gap_fraction"]) for row in per_sequence],
            "Sequence",
            "Gap fraction",
            "MSATK Missingness by Sequence",
        ),
    )
    save(
        "gap_fraction_by_site.png",
        lambda: _line(
            plt,
            [int(row["site_index"]) for row in per_site],
            [float(row["gap_fraction"]) for row in per_site],
            "Site",
            "Gap fraction",
            "MSATK Gap Fraction by Site",
        ),
    )
    save(
        "entropy_by_site.png",
        lambda: _line(
            plt,
            [int(row["site_index"]) for row in per_site],
            [float(row["entropy"]) for row in per_site],
            "Site",
            "Entropy",
            "MSATK Entropy by Site",
        ),
    )
    save(
        "conservation_by_site.png",
        lambda: _line(
            plt,
            [int(row["site_index"]) for row in per_site],
            [float(row["conservation_score"]) for row in per_site],
            "Site",
            "Conservation score",
            "MSATK Conservation by Site",
        ),
    )
    save(
        "composition_summary.png",
        lambda: _bar(
            plt,
            [
                str(row.get("character", row.get("amino_acid", row.get("codon", ""))))
                for row in composition
            ],
            [float(row.get("fraction", row.get("frequency", 0.0))) for row in composition],
            "Character",
            "Fraction",
            "MSATK Composition Summary",
        ),
    )
    save(
        "pairwise_identity_heatmap.png",
        lambda: _heatmap(plt, identities, "MSATK Pairwise Identity Heatmap"),
    )
    distances = [
        1.0 - identities[i][j]
        for i in range(len(identities))
        for j in range(i + 1, len(identities))
    ]
    save(
        "pairwise_distance_distribution.png",
        lambda: _hist(plt, distances, "Distance", "Count", "MSATK Pairwise Distance Distribution"),
    )
    pairwise_identities = [
        identities[i][j] for i in range(len(identities)) for j in range(i + 1, len(identities))
    ]
    save(
        "pairwise_identity_distribution.png",
        lambda: _hist(
            plt, pairwise_identities, "Identity", "Count", "MSATK Pairwise Identity Distribution"
        ),
    )
    if embeddings:
        save(
            "sequence_embedding.png",
            lambda: _scatter(
                plt,
                [float(row["embedding_1"]) for row in embeddings],
                [float(row["embedding_2"]) for row in embeddings],
                "Embedding 1",
                "Embedding 2",
                "MSATK Sequence Embedding",
            ),
        )
    return generated


def plot_single_series(
    path: str | Path, x: list[int], y: list[float], xlabel: str, ylabel: str, title: str
) -> bool:
    """Write a single line plot. Returns False when matplotlib is unavailable."""

    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 4.8), dpi=150)
    _line(plt, x, y, xlabel, ylabel, title)
    plt.tight_layout()
    plt.savefig(target)
    plt.close()
    return True


def plot_single_bar(
    path: str | Path, labels: list[str], values: list[float], xlabel: str, ylabel: str, title: str
) -> bool:
    """Write a single bar plot. Returns False when matplotlib is unavailable."""

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


def _bar(plt, labels: list[str], values: list[float], xlabel: str, ylabel: str, title: str) -> None:
    plt.bar(range(len(values)), values, color="#0f766e")
    if len(labels) <= 30:
        plt.xticks(range(len(labels)), labels, rotation=60, ha="right", fontsize=7)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def _line(plt, x: list[int], y: list[float], xlabel: str, ylabel: str, title: str) -> None:
    plt.plot(x, y, color="#0f766e", linewidth=1.6)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def _heatmap(plt, matrix: list[list[float]], title: str) -> None:
    plt.imshow(matrix, cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(label="Identity")
    plt.title(title)


def _hist(plt, values: list[float], xlabel: str, ylabel: str, title: str) -> None:
    plt.hist(values or [0.0], bins=min(30, max(5, len(values))), color="#0f766e", edgecolor="white")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def _scatter(plt, x: list[float], y: list[float], xlabel: str, ylabel: str, title: str) -> None:
    plt.scatter(x, y, color="#0f766e", s=36)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
