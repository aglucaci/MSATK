"""Numerical tensor encodings of multiple sequence alignments.

All encoders return ``numpy`` arrays with sequences on axis 0 and aligned sites
on axis 1, plus a list describing the channel semantics. Row and column order is
exactly the input order unless a transform (see :mod:`msatk.encoding.transforms`)
was applied, so ``tensor[i, j]`` always maps back to ``(row i, site j+1)``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from msatk.constants import HYDROPATHY_KD
from msatk.core.alphabet import is_nucleotide
from msatk.core.encoded import (
    STATE_AMBIGUOUS,
    STATE_GAP,
    STATE_MISSING,
    EncodedAlignment,
    encode_alignment,
)
from msatk.encoding.vocab import Vocabulary, ambiguity_expansion, full_vocabulary
from msatk.models import Alignment, SequenceRecord
from msatk.stats.site import site_arrays

ENCODINGS = ("categorical", "onehot", "physicochemical", "features", "rgb")
AMBIGUITY_POLICIES = ("missing", "fractional", "separate")


@dataclass
class EncodedTensor:
    """An encoded alignment tensor plus everything needed to interpret it."""

    array: np.ndarray
    encoding: str
    channels: list[str]
    ids: list[str]
    molecule_type: str
    lossless: bool
    details: dict[str, object] = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(int(s) for s in self.array.shape)


def _enc(alignment: Alignment | EncodedAlignment, molecule_type: str) -> EncodedAlignment:
    return (
        alignment
        if isinstance(alignment, EncodedAlignment)
        else encode_alignment(alignment, molecule_type)
    )


def categorical(
    alignment: Alignment, molecule_type: str, vocabulary: Vocabulary | None = None
) -> EncodedTensor:
    """Integer label per cell over a fixed full vocabulary (exactly reversible).

    Symbols are case-folded. Characters outside the vocabulary map to
    ``vocabulary.other_index`` and make the encoding non-lossless.
    """

    vocab = vocabulary or full_vocabulary(molecule_type)
    lut = np.full(256, vocab.other_index, dtype=np.uint8)
    for symbol, idx in vocab.index().items():
        lut[ord(symbol)] = idx
        lut[ord(symbol.lower())] = idx
    seqs = alignment.padded_sequences()
    raw = np.frombuffer("".join(seqs).encode("ascii", errors="replace"), dtype=np.uint8)
    labels = lut[raw].reshape(alignment.n_sequences, alignment.length)
    other = int((labels == vocab.other_index).sum())
    return EncodedTensor(
        labels,
        "categorical",
        [*vocab.symbols, "<other>"],
        alignment.ids,
        molecule_type,
        lossless=other == 0 and alignment.is_rectangular,
        details={"vocabulary": vocab.describe(), "other_cells": other, "dtype": "uint8"},
    )


def decode_categorical(
    labels: np.ndarray, ids: list[str], vocabulary: Vocabulary, other: str = "?"
) -> Alignment:
    """Invert :func:`categorical` (upper-case symbols)."""

    table = np.array([*vocabulary.symbols, other])
    rows = ["".join(table[row]) for row in labels]
    return Alignment(tuple(SequenceRecord(i, s) for i, s in zip(ids, rows)))


def onehot(
    alignment: Alignment | EncodedAlignment,
    molecule_type: str,
    ambiguity: str = "missing",
    dtype: str = "float32",
) -> EncodedTensor:
    """One-hot (or soft) canonical-residue channels plus explicit gap / missing channels.

    ``ambiguity``:
      * ``missing``    – ambiguity codes are encoded in the ``<missing>`` channel.
      * ``fractional`` – ambiguity codes (and N/X) are spread uniformly over the
        compatible canonical residues (e.g. ``R`` → A=0.5, G=0.5); ``?`` stays missing.
      * ``separate``   – ambiguity codes get their own ``<ambiguous>`` channel.

    Every non-padding cell sums to exactly 1 across channels.
    """

    if ambiguity not in AMBIGUITY_POLICIES:
        raise ValueError(f"ambiguity must be one of {AMBIGUITY_POLICIES}")
    enc = _enc(alignment, molecule_type)
    k = enc.k
    channels = list(enc.policy.canonical) + ["<gap>"]
    if ambiguity == "separate":
        channels.append("<ambiguous>")
    channels.append("<missing>")
    x = np.zeros((enc.n, enc.width, len(channels)), dtype=np.float32)
    for c in range(k):
        x[:, :, c] = enc.residue_index == c
    x[:, :, k] = enc.gap_mask
    amb = enc.state == STATE_AMBIGUOUS
    missing = (enc.state == STATE_MISSING) | (enc.state >= 4)  # missing, stop, invalid
    if ambiguity == "separate":
        x[:, :, k + 1] = amb
        x[:, :, k + 2] = missing
    elif ambiguity == "missing":
        x[:, :, k + 1] = amb | missing
    else:
        expansion = ambiguity_expansion(enc.policy.molecule_type)
        canon = {ch: i for i, ch in enumerate(enc.policy.canonical)}
        handled = np.zeros((enc.n, enc.width), dtype=bool)
        candidates = ~enc.residue_mask & ~enc.gap_mask
        for code, bases in expansion.items():
            cells = candidates & (enc.chars == ord(code))
            if not cells.any():
                continue
            for b in bases:
                x[:, :, canon[b]][cells] = 1.0 / len(bases)
            handled |= cells
        x[:, :, k + 1] = (amb | missing) & ~handled
    out = x.astype(dtype)
    return EncodedTensor(
        out,
        "onehot",
        channels,
        enc.alignment.ids,
        enc.policy.molecule_type,
        lossless=False,
        details={
            "ambiguity_policy": ambiguity,
            "note": "Lossless for canonical residues and gaps; ambiguity codes collapse per policy.",
            "dtype": dtype,
        },
    )


_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
_VOLUME = {  # Zamyatnin (1972) residue volumes, Å^3
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4,
    "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9,
    "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}  # fmt: skip
_POLAR = set("STNQCYHKRDE")
_AROMATIC = set("FWYH")


def physicochemical(
    alignment: Alignment | EncodedAlignment, molecule_type: str = "protein"
) -> EncodedTensor:
    """Biochemically meaningful continuous channels for protein alignments.

    Channels: Kyte–Doolittle hydropathy (/4.5), net charge at pH ~7, residue volume
    (min–max scaled), polar flag, aromatic flag, gap mask, missing mask. Non-residue
    cells are zero in property channels and flagged in the masks.
    """

    enc = _enc(alignment, molecule_type)
    if is_nucleotide(enc.policy.molecule_type):
        raise ValueError("Physicochemical channels are defined for protein alignments only.")
    alphabet = enc.policy.canonical
    vmin, vmax = min(_VOLUME.values()), max(_VOLUME.values())
    props = np.array(
        [
            [
                HYDROPATHY_KD[a] / 4.5,
                _CHARGE.get(a, 0.0),
                (_VOLUME[a] - vmin) / (vmax - vmin),
                float(a in _POLAR),
                float(a in _AROMATIC),
            ]
            for a in alphabet
        ],
        dtype=np.float32,
    )
    x = np.zeros((enc.n, enc.width, 7), dtype=np.float32)
    mask = enc.residue_mask
    x[:, :, :5][mask] = props[enc.residue_index[mask]]
    x[:, :, 5] = enc.gap_mask
    x[:, :, 6] = ~mask & ~enc.gap_mask
    return EncodedTensor(
        x,
        "physicochemical",
        ["hydropathy_kd", "charge", "volume", "polar", "aromatic", "<gap>", "<non_residue>"],
        enc.alignment.ids,
        enc.policy.molecule_type,
        lossless=False,
        details={"scales": "KD/4.5; charge K,R=+1 D,E=-1 H=+0.1; volume Zamyatnin 1972 min-max"},
    )


def feature_channels(alignment: Alignment | EncodedAlignment, molecule_type: str) -> EncodedTensor:
    """Per-cell QC/feature channels combining cell state and column statistics."""

    enc = _enc(alignment, molecule_type)
    sites = site_arrays(enc)
    n, w = enc.n, enc.width
    consensus = enc.consensus_index
    x = np.zeros((n, w, 7), dtype=np.float32)
    x[:, :, 0] = enc.residue_mask
    x[:, :, 1] = enc.state == STATE_GAP
    x[:, :, 2] = (enc.state == STATE_MISSING) | (enc.state >= 4)
    x[:, :, 3] = enc.state == STATE_AMBIGUOUS
    x[:, :, 4] = (enc.residue_index == consensus[None, :]) & enc.residue_mask
    x[:, :, 5] = np.nan_to_num(sites["normalized_entropy"], nan=0.0)[None, :]
    x[:, :, 6] = (sites["observed"] / max(n, 1))[None, :]
    return EncodedTensor(
        x,
        "features",
        [
            "residue",
            "gap",
            "missing",
            "ambiguous",
            "matches_consensus",
            "site_normalized_entropy",
            "site_occupancy",
        ],
        enc.alignment.ids,
        enc.policy.molecule_type,
        lossless=False,
        details={"note": "Site channels are broadcast along the sequence axis."},
    )


def encode(
    alignment: Alignment, molecule_type: str, encoding: str = "onehot", **kwargs: object
) -> EncodedTensor:
    if encoding == "categorical":
        return categorical(alignment, molecule_type)
    if encoding == "onehot":
        return onehot(alignment, molecule_type, **kwargs)  # type: ignore[arg-type]
    if encoding == "physicochemical":
        return physicochemical(alignment, molecule_type)
    if encoding == "features":
        return feature_channels(alignment, molecule_type)
    if encoding == "rgb":
        from msatk.encoding.images import rgb_tensor

        return rgb_tensor(alignment, molecule_type, **kwargs)  # type: ignore[arg-type]
    raise ValueError(f"Unknown encoding {encoding!r}; choose from {ENCODINGS}.")
