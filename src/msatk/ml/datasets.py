"""Framework-agnostic dataset wrappers for exported MSATK tensors.

:class:`AlignmentTensorDataset` implements ``__len__``/``__getitem__`` and returns
NumPy arrays, so it can be passed directly to ``torch.utils.data.DataLoader`` or
wrapped with ``tf.data.Dataset.from_generator``. Framework adapters import the
framework lazily and raise a clear error when it is not installed.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from msatk.encoding.export import load_tensor
from msatk.exceptions import MissingDependencyError


class AlignmentTensorDataset:
    """Items from one or more MSATK exports.

    * ``tiles``/``pad`` exports: one item per tile, ``{"x", "mask", "tile", "source"}``.
      Several files can be combined when their tiles have the same shape and
      channels, e.g. alignments of different sizes exported with the same
      ``--height``/``--width``.
    * ``native`` exports: one item per file, ``{"x", "mask", "source"}``, at the
      alignment's own size. Batch them with :func:`pad_batch`.

    ``source`` is the index of the file in ``npz_paths``.
    """

    def __init__(
        self, npz_paths: str | Path | Sequence[str | Path], channels_first: bool = False
    ) -> None:
        paths = [npz_paths] if isinstance(npz_paths, (str, Path)) else list(npz_paths)
        if not paths:
            raise ValueError("No exports given.")
        self.paths = [Path(p) for p in paths]
        self.channels_first = channels_first
        loaded = [load_tensor(p) for p in self.paths]
        self.manifests = [m for _, m in loaded]
        self.manifest = self.manifests[0]
        kinds = {"tiles" in a for a, _ in loaded}
        if len(kinds) > 1:
            raise ValueError("Cannot mix native exports with tiles/pad exports in one dataset.")
        self.windowed = kinds.pop()
        channels = {json.dumps(m.get("channels")) for m in self.manifests}
        if len(channels) > 1:
            raise ValueError("Exports use different encodings/channels and cannot be combined.")
        self._items: list[tuple[int, int]] = []
        self._arrays = [a for a, _ in loaded]
        if self.windowed:
            shapes = {a["x"].shape[1:] for a in self._arrays}
            if len(shapes) > 1:
                raise ValueError(
                    f"Tile shapes differ between exports ({sorted(shapes)}); export them with "
                    "the same --height and --width."
                )
            for k, a in enumerate(self._arrays):
                self._items += [(k, t) for t in range(a["x"].shape[0])]
        else:
            self._items = [(k, 0) for k in range(len(self._arrays))]
        # Backwards-compatible attributes for single-file use.
        first = self._arrays[0]
        self.arrays = first
        self.x = first["x"] if self.windowed else first["x"][None, ...]
        self.mask = first.get("mask")

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, index: int) -> dict[str, Any]:
        k, t = self._items[index]
        arrays = self._arrays[k]
        if self.windowed:
            x = arrays["x"][t]
            mask = arrays["mask"][t]
        else:
            x = arrays["x"]
            mask = np.ones(x.shape[:2], dtype=bool)
        if self.channels_first and x.ndim == 3:
            x = np.moveaxis(x, -1, 0)
        item: dict[str, Any] = {"x": x, "mask": mask, "source": k}
        if self.windowed:
            item["tile"] = arrays["tiles"][t]
        return item


def pad_batch(
    items: Sequence[dict[str, Any]], multiple_of: int = 1, channels_first: bool | None = None
) -> dict[str, Any]:
    """Collate items of different sizes by padding to the largest in the batch.

    For models that accept variable input sizes (e.g. CNNs with adaptive pooling).
    Returns NumPy arrays ``x`` ``(B, rows, sites[, C])`` (or ``(B, C, rows, sites)``
    when the items are channels-first), ``mask`` ``(B, rows, sites)`` with False on
    padding, and ``source``. ``multiple_of`` rounds the padded size up, e.g. to a
    patch or downsampling factor. Use as ``DataLoader(..., collate_fn=pad_batch)``
    and convert with :func:`to_torch`.
    """

    if not items:
        raise ValueError("Empty batch.")
    masks = [np.asarray(it["mask"], dtype=bool) for it in items]
    xs = [np.asarray(it["x"]) for it in items]
    if channels_first is None:
        channels_first = xs[0].ndim == 3 and xs[0].shape[1:] == masks[0].shape
    m = max(1, multiple_of)

    def up(v: int) -> int:
        return -(-v // m) * m

    rows = up(max(mk.shape[0] for mk in masks))
    sites = up(max(mk.shape[1] for mk in masks))
    batch_mask = np.zeros((len(items), rows, sites), dtype=bool)
    if xs[0].ndim == 2:
        out = np.zeros((len(items), rows, sites), dtype=xs[0].dtype)
    elif channels_first:
        out = np.zeros((len(items), xs[0].shape[0], rows, sites), dtype=xs[0].dtype)
    else:
        out = np.zeros((len(items), rows, sites, xs[0].shape[-1]), dtype=xs[0].dtype)
    for b, (x, mk) in enumerate(zip(xs, masks)):
        r, c = mk.shape
        if xs[0].ndim == 3 and channels_first:
            out[b, :, :r, :c] = x
        else:
            out[b, :r, :c] = x
        batch_mask[b, :r, :c] = mk
    return {"x": out, "mask": batch_mask, "source": np.array([it.get("source", 0) for it in items])}


def to_torch(array: np.ndarray) -> Any:
    """Convert a NumPy array to a ``torch.Tensor`` (requires PyTorch)."""

    try:
        import torch
    except ImportError as exc:
        raise MissingDependencyError("PyTorch is not installed: pip install torch") from exc
    return torch.from_numpy(np.ascontiguousarray(array))


def to_tensorflow(array: np.ndarray) -> Any:
    """Convert a NumPy array to a ``tf.Tensor`` (requires TensorFlow)."""

    try:
        import tensorflow as tf  # type: ignore
    except ImportError as exc:
        raise MissingDependencyError("TensorFlow is not installed: pip install tensorflow") from exc
    return tf.convert_to_tensor(array)
