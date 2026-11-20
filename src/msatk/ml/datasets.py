"""Framework-agnostic dataset wrappers for exported MSATK tensors.

:class:`AlignmentTensorDataset` implements ``__len__``/``__getitem__`` and returns
NumPy arrays, so it can be passed directly to ``torch.utils.data.DataLoader`` or
wrapped with ``tf.data.Dataset.from_generator``. Framework adapters import the
framework lazily and raise a clear error when it is not installed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from msatk.encoding.export import load_tensor
from msatk.exceptions import MissingDependencyError


class AlignmentTensorDataset:
    """Iterate over the tiles of a windowed export (or the single tensor of an unwindowed one)."""

    def __init__(self, npz_path: str | Path, channels_first: bool = False) -> None:
        self.arrays, self.manifest = load_tensor(npz_path)
        x = self.arrays["x"]
        self.windowed = "tiles" in self.arrays
        self.x = x if self.windowed else x[None, ...]
        self.mask = self.arrays.get("mask")
        self.channels_first = channels_first

    def __len__(self) -> int:
        return int(self.x.shape[0])

    def __getitem__(self, index: int) -> dict[str, Any]:
        x = self.x[index]
        if self.channels_first and x.ndim == 3:
            x = np.moveaxis(x, -1, 0)
        item: dict[str, Any] = {"x": x}
        if self.mask is not None:
            item["mask"] = self.mask[index]
        if self.windowed:
            item["tile"] = self.arrays["tiles"][index]
        return item


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
