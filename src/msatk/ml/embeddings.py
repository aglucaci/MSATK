"""Low-dimensional sequence embeddings (PCA, UMAP, t-SNE).

Embeddings are optional and never computed on the default ``profile`` path.
The algorithm that actually ran is always reported: MSATK raises
:class:`~msatk.exceptions.MissingDependencyError` instead of silently falling
back to a different method.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from msatk.core.detect import detect_molecule_type
from msatk.core.encoded import encode_alignment
from msatk.exceptions import MissingDependencyError
from msatk.models import Alignment

EMBEDDING_METHODS = ("pca", "umap", "tsne")


def one_hot_features(
    alignment: Alignment, molecule_type: str | None = None
) -> tuple[list[str], np.ndarray]:
    """Flattened one-hot canonical-residue features plus a gap channel, shape ``(N, L*(K+1))``."""

    enc = encode_alignment(alignment, molecule_type or detect_molecule_type(alignment))
    channels = list(enc.policy.canonical) + ["-"]
    feats = np.zeros((enc.n, enc.width, len(channels)), dtype=np.float32)
    for c in range(enc.k):
        feats[:, :, c] = enc.residue_index == c
    feats[:, :, -1] = enc.gap_mask
    names = [f"site_{i + 1}_{ch}" for i in range(enc.width) for ch in channels]
    return names, feats.reshape(enc.n, -1)


def kmer_features(alignment: Alignment, k: int = 3) -> tuple[list[str], np.ndarray]:
    """Alignment-free k-mer frequency features computed on ungapped sequences."""

    ungapped = ["".join(ch for ch in seq.upper() if ch not in "-.") for seq in alignment.sequences]
    kmers = sorted({s[i : i + k] for s in ungapped for i in range(max(len(s) - k + 1, 0))})
    index = {kmer: j for j, kmer in enumerate(kmers)}
    out = np.zeros((len(ungapped), len(kmers)), dtype=np.float32)
    for row, s in enumerate(ungapped):
        counts = Counter(s[i : i + k] for i in range(max(len(s) - k + 1, 0)))
        total = sum(counts.values())
        for kmer, count in counts.items():
            out[row, index[kmer]] = count / total if total else 0.0
    return [f"kmer_{kmer}" for kmer in kmers], out


def sequence_embeddings(
    alignment: Alignment,
    method: str = "pca",
    representation: str = "onehot",
    molecule_type: str | None = None,
    random_state: int = 7,
) -> list[dict[str, object]]:
    """Return a two-dimensional embedding of each sequence."""

    method = method.lower()
    if method not in EMBEDDING_METHODS:
        raise ValueError(f"Unknown embedding method {method!r}; choose from {EMBEDDING_METHODS}.")
    _, x = (
        kmer_features(alignment)
        if representation == "kmer"
        else one_hot_features(alignment, molecule_type)
    )
    n = x.shape[0]
    explained: list[float | None] = [None, None]
    if n < 3 or x.shape[1] == 0:
        coords = np.zeros((n, 2))
        method_used = "none (fewer than 3 sequences)"
    elif method == "pca":
        coords, explained = _pca(x)
        method_used = "pca"
    else:
        try:
            if method == "umap":
                import umap  # type: ignore

                reducer = umap.UMAP(
                    n_components=2, random_state=random_state, n_neighbors=min(15, n - 1)
                )
            else:
                from sklearn.manifold import TSNE

                reducer = TSNE(
                    n_components=2,
                    init="random",
                    perplexity=float(max(1, min(30, n - 1))),
                    learning_rate="auto",
                    random_state=random_state,
                )
        except ImportError as exc:
            extra = "umap-learn" if method == "umap" else "scikit-learn"
            raise MissingDependencyError(
                f"Embedding method {method!r} requires {extra}. Install with: "
                f"pip install 'msatk[embed]' — or use --embeddings pca."
            ) from exc
        coords = np.asarray(reducer.fit_transform(x))
        method_used = method
    return [
        {
            "sequence_id": seq_id,
            "embedding_1": float(xy[0]),
            "embedding_2": float(xy[1]),
            "method": method_used,
            "representation": representation,
            "explained_variance_ratio_1": explained[0],
            "explained_variance_ratio_2": explained[1],
        }
        for seq_id, xy in zip(alignment.ids, coords)
    ]


def _pca(x: np.ndarray) -> tuple[np.ndarray, list[float | None]]:
    """Exact PCA via SVD (numpy only; deterministic sign convention)."""

    centered = x - x.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    signs = np.sign(vt[np.arange(vt.shape[0]), np.abs(vt).argmax(axis=1)])
    signs[signs == 0] = 1
    coords = (u * s) * signs
    var = s**2
    total = var.sum()
    coords2 = np.zeros((x.shape[0], 2))
    k = min(2, coords.shape[1])
    coords2[:, :k] = coords[:, :k]
    ratios: list[float | None] = [
        float(var[i] / total) if total and i < var.size else None for i in range(2)
    ]
    return coords2, ratios
