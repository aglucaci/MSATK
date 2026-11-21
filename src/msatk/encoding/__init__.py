"""Machine-learning representations of multiple sequence alignments.

* :mod:`~msatk.encoding.tensors` – categorical, one-hot, physicochemical and feature tensors.
* :mod:`~msatk.encoding.images` – RGB residue images and lossless label images.
* :mod:`~msatk.encoding.transforms` – coordinate-preserving crops, reordering and windowing.
* :mod:`~msatk.encoding.ordering` – recorded row orderings (reference, clustering, user lists).
* :mod:`~msatk.encoding.export` – NPZ/PNG writers with self-describing manifests.
"""

from .export import export_image, export_tensor, load_tensor
from .images import label_image, rgb_image, rgb_tensor
from .ordering import ORDER_METHODS, row_order
from .tensors import (
    ENCODINGS,
    EncodedTensor,
    categorical,
    decode_categorical,
    encode,
    feature_channels,
    onehot,
    physicochemical,
)
from .transforms import CoordinateMap, Windows, sliding_windows, subset_alignment
from .vocab import NUCLEOTIDE_FULL, PROTEIN_FULL, Vocabulary, full_vocabulary

__all__ = [
    "ENCODINGS",
    "ORDER_METHODS",
    "CoordinateMap",
    "EncodedTensor",
    "NUCLEOTIDE_FULL",
    "PROTEIN_FULL",
    "Vocabulary",
    "Windows",
    "categorical",
    "decode_categorical",
    "encode",
    "export_image",
    "export_tensor",
    "feature_channels",
    "full_vocabulary",
    "label_image",
    "load_tensor",
    "onehot",
    "physicochemical",
    "rgb_image",
    "rgb_tensor",
    "row_order",
    "sliding_windows",
    "subset_alignment",
]
