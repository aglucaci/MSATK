"""Alignment-to-image rendering for visual inspection and image-model pipelines.

Two families of images are produced and must not be confused:

* **RGB residue images** (:func:`rgb_image`) – colours are arbitrary labels for
  humans or for RGB-pretrained models. Integer upscaling (each cell → an
  ``s × s`` block) is lossless; resizing to a fixed size is *lossy* and recorded.
* **Label images** (:func:`label_image`) – single-channel uint8 images whose pixel
  value is the categorical vocabulary index, i.e. an exactly reversible encoding.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from msatk.encoding.tensors import EncodedTensor, categorical
from msatk.models import Alignment
from msatk.viz.themes import OTHER_COLOR, hex_to_rgb, load_theme


def rgb_image(
    alignment: Alignment,
    molecule_type: str,
    theme: str | dict[str, str] | Path | None = "auto",
    scale: int = 1,
) -> np.ndarray:
    """``(N*scale, L*scale, 3)`` uint8 RGB image of the alignment."""

    mapping = load_theme(theme, molecule_type)
    other = hex_to_rgb(OTHER_COLOR)
    lut = np.tile(np.array(other, dtype=np.uint8), (256, 1))
    if not [k for k in mapping if k not in "-.?"]:
        # Grayscale theme: residues dark, gaps light.
        lut[:] = (60, 60, 60)
    for symbol, color in mapping.items():
        if len(symbol) != 1:
            continue
        rgb = np.array(hex_to_rgb(color), dtype=np.uint8)
        lut[ord(symbol.upper())] = rgb
        lut[ord(symbol.lower())] = rgb
    raw = np.frombuffer(
        "".join(alignment.padded_sequences()).encode("ascii", errors="replace"), dtype=np.uint8
    ).reshape(alignment.n_sequences, alignment.length)
    img = lut[raw]
    if scale > 1:
        img = np.repeat(np.repeat(img, scale, axis=0), scale, axis=1)
    return img


def rgb_tensor(
    alignment: Alignment,
    molecule_type: str,
    theme: str | dict[str, str] | Path | None = "auto",
) -> EncodedTensor:
    """``(N, L, 3)`` uint8 RGB tensor (one pixel per cell) for RGB-pretrained models.

    The colours are display labels, not measurements: two residues with similar
    colours are not biologically similar, and residues sharing a colour in the
    theme become indistinguishable. Use ``onehot`` for scientific work.
    """

    img = rgb_image(alignment, molecule_type, theme=theme, scale=1)
    return EncodedTensor(
        array=img,
        encoding="rgb",
        channels=["R", "G", "B"],
        ids=alignment.ids,
        molecule_type=molecule_type,
        lossless=False,
        details={
            "theme": theme if isinstance(theme, str) or theme is None else "custom",
            "pixels_per_cell": 1,
            "note": "RGB colours are arbitrary display labels; symbols that share a colour "
            "cannot be told apart. Provided for compatibility with RGB-pretrained models.",
        },
    )


def label_image(alignment: Alignment, molecule_type: str) -> tuple[np.ndarray, list[str]]:
    """Single-channel uint8 image of vocabulary indices (lossless) and the index legend."""

    tensor = categorical(alignment, molecule_type)
    return tensor.array.astype(np.uint8), tensor.channels


def resize_nearest(img: np.ndarray, height: int, width: int) -> np.ndarray:
    """Nearest-neighbour resize (no interpolation of colours). Lossy when shrinking."""

    rows = (np.arange(height) * img.shape[0] / height).astype(int)
    cols = (np.arange(width) * img.shape[1] / width).astype(int)
    return img[rows][:, cols]


def save_image(img: np.ndarray, path: str | Path) -> Path:
    """Write a uint8 array as PNG/TIFF (via matplotlib/Pillow) or NPY."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.suffix.lower() == ".npy":
        np.save(target, img)
        return target
    try:
        from PIL import Image
    except ImportError:  # pragma: no cover - Pillow ships with matplotlib
        import matplotlib.pyplot as plt

        plt.imsave(target, img, cmap="gray" if img.ndim == 2 else None, vmin=0, vmax=255)
        return target
    Image.fromarray(img).save(target)
    return target
