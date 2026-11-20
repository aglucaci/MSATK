"""Vectorised alignment representation shared by statistics and ML encoders."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from msatk.core.alphabet import AlphabetPolicy, policy_for
from msatk.models import Alignment

STATE_RESIDUE = 0
STATE_AMBIGUOUS = 1
STATE_MISSING = 2
STATE_GAP = 3
STATE_STOP = 4
STATE_INVALID = 5
STATE_NAMES = ("residue", "ambiguous", "missing", "gap", "stop", "invalid")


@dataclass
class EncodedAlignment:
    """Character, state and residue-index matrices for an alignment.

    Attributes
    ----------
    chars:
        ``(N, L)`` uint8 matrix of normalised (upper-case, synonym-mapped) ASCII codes.
    state:
        ``(N, L)`` uint8 matrix of cell states (see :mod:`msatk.core.alphabet`).
    residue_index:
        ``(N, L)`` int16 matrix with the canonical-alphabet index of each residue, or ``-1``.
    """

    alignment: Alignment
    policy: AlphabetPolicy
    chars: np.ndarray
    state: np.ndarray
    residue_index: np.ndarray
    original_lengths: np.ndarray
    _cache: dict[str, object] = field(default_factory=dict, repr=False)

    @property
    def n(self) -> int:
        return int(self.chars.shape[0])

    @property
    def width(self) -> int:
        return int(self.chars.shape[1])

    @property
    def k(self) -> int:
        return self.policy.alphabet_size

    @cached_property
    def residue_mask(self) -> np.ndarray:
        return self.state == STATE_RESIDUE

    @cached_property
    def gap_mask(self) -> np.ndarray:
        return self.state == STATE_GAP

    @cached_property
    def pad_mask(self) -> np.ndarray:
        cols = np.arange(self.width)[None, :]
        return cols >= self.original_lengths[:, None]

    @cached_property
    def site_counts(self) -> np.ndarray:
        """``(L, K)`` counts of each canonical residue per column."""

        counts = np.zeros((self.width, self.k), dtype=np.int64)
        for index in range(self.k):
            counts[:, index] = (self.residue_index == index).sum(axis=0)
        return counts

    @cached_property
    def sequence_counts(self) -> np.ndarray:
        """``(N, K)`` counts of each canonical residue per sequence."""

        counts = np.zeros((self.n, self.k), dtype=np.int64)
        for index in range(self.k):
            counts[:, index] = (self.residue_index == index).sum(axis=1)
        return counts

    @cached_property
    def consensus_index(self) -> np.ndarray:
        """Majority canonical residue per column (``-1`` for unobserved columns).

        Ties are broken by canonical alphabet order, which keeps results deterministic.
        """

        counts = self.site_counts
        idx = counts.argmax(axis=1).astype(np.int16)
        idx[counts.sum(axis=1) == 0] = -1
        return idx

    def consensus_string(self, unobserved: str = "-") -> str:
        alphabet = self.policy.canonical
        return "".join(alphabet[i] if i >= 0 else unobserved for i in self.consensus_index)

    def row_strings(self) -> list[str]:
        return [bytes(row).decode("ascii") for row in self.chars]


def encode_alignment(alignment: Alignment, molecule_type: str) -> EncodedAlignment:
    """Vectorise an alignment under the character policy for ``molecule_type``."""

    policy = policy_for(molecule_type)
    seqs = alignment.padded_sequences()
    width = alignment.length
    raw = "".join(seqs).encode("ascii", errors="replace")
    chars = np.frombuffer(raw, dtype=np.uint8).reshape(len(seqs), width).copy()

    # Build 256-entry lookup tables once.
    norm = np.arange(256, dtype=np.uint8)
    state_lut = np.full(256, STATE_INVALID, dtype=np.uint8)
    index_lut = np.full(256, -1, dtype=np.int16)
    canonical = {ch: i for i, ch in enumerate(policy.canonical)}
    state_codes = {
        "residue": STATE_RESIDUE,
        "ambiguous": STATE_AMBIGUOUS,
        "missing": STATE_MISSING,
        "gap": STATE_GAP,
        "stop": STATE_STOP,
        "invalid": STATE_INVALID,
    }
    for code in range(32, 127):
        ch = chr(code)
        normalized = ch if ch in policy.gaps else policy.normalize(ch)
        norm[code] = ord(normalized) if len(normalized) == 1 else code
        state_lut[code] = state_codes[policy.state(ch)]
        if state_lut[code] == STATE_RESIDUE:
            index_lut[code] = canonical[normalized]
    state = state_lut[chars]
    residue_index = index_lut[chars]
    chars = norm[chars]
    lengths = np.array([len(record.sequence) for record in alignment.records], dtype=np.int64)
    return EncodedAlignment(alignment, policy, chars, state, residue_index, lengths)
