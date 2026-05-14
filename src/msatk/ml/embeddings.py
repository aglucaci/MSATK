"""Sequence embedding helpers."""

from __future__ import annotations

from collections import Counter

from msatk.constants import AA_CANONICAL, DNA_BASES
from msatk.models import Alignment


def one_hot_features(alignment: Alignment) -> tuple[list[str], list[list[float]]]:
    alphabet = sorted(
        (DNA_BASES | AA_CANONICAL | {"-"})
        & {ch.upper() for seq in alignment.sequences for ch in seq}
    )
    if not alphabet:
        alphabet = ["-"]
    columns = [(index, char) for index in range(alignment.length) for char in alphabet]
    rows: list[list[float]] = []
    for seq in alignment.padded_sequences():
        rows.append([1.0 if seq[index].upper() == char else 0.0 for index, char in columns])
    return [f"site_{idx + 1}_{char}" for idx, char in columns], rows


def kmer_features(alignment: Alignment, k: int = 3) -> tuple[list[str], list[list[float]]]:
    kmers = sorted(
        {
            seq[i : i + k].upper()
            for seq in alignment.sequences
            for i in range(max(len(seq) - k + 1, 0))
            if "-" not in seq[i : i + k]
        }
    )
    rows: list[list[float]] = []
    for seq in alignment.sequences:
        counts = Counter(seq[i : i + k].upper() for i in range(max(len(seq) - k + 1, 0)))
        total = sum(counts.values())
        rows.append([counts.get(kmer, 0) / total if total else 0.0 for kmer in kmers])
    return [f"kmer_{kmer}" for kmer in kmers], rows


def sequence_embeddings(
    alignment: Alignment, method: str = "pca", representation: str = "onehot"
) -> list[dict[str, object]]:
    """Return sequence embeddings. Uses scikit-learn when available, otherwise first two features."""

    _, features = (
        kmer_features(alignment) if representation == "kmer" else one_hot_features(alignment)
    )
    method = method.lower()
    coords: list[list[float]]
    try:
        import numpy as np
        from sklearn.decomposition import PCA
        from sklearn.manifold import TSNE

        x = np.asarray(features, dtype=float)
        if x.shape[0] == 1:
            coords = [[0.0, 0.0]]
        elif method == "tsne":
            perplexity = max(1, min(30, x.shape[0] - 1))
            coords = (
                TSNE(n_components=2, init="random", perplexity=perplexity, learning_rate="auto")
                .fit_transform(x)
                .tolist()
            )
        elif method == "umap":
            try:
                import umap  # type: ignore

                coords = umap.UMAP(n_components=2, random_state=7).fit_transform(x).tolist()
            except Exception:
                coords = PCA(n_components=2).fit_transform(x).tolist()
        else:
            coords = PCA(n_components=2).fit_transform(x).tolist()
    except Exception:
        coords = [[row[0] if row else 0.0, row[1] if len(row) > 1 else 0.0] for row in features]
    return [
        {
            "sequence_id": seq_id,
            "embedding_1": xy[0],
            "embedding_2": xy[1],
            "method": method,
            "representation": representation,
        }
        for seq_id, xy in zip(alignment.ids, coords)
    ]
