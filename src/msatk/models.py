"""Core data structures used by MSATK."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SequenceRecord:
    """One aligned sequence."""

    id: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class Alignment:
    """A multiple sequence alignment."""

    records: tuple[SequenceRecord, ...]
    source: str = ""
    fmt: str = "unknown"

    def __post_init__(self) -> None:
        if not self.records:
            raise ValueError("MSATK requires at least one sequence.")

    @property
    def ids(self) -> list[str]:
        return [record.id for record in self.records]

    @property
    def sequences(self) -> list[str]:
        return [record.sequence for record in self.records]

    @property
    def length(self) -> int:
        return max(len(record.sequence) for record in self.records)

    @property
    def n_sequences(self) -> int:
        return len(self.records)

    @property
    def is_rectangular(self) -> bool:
        lengths = {len(record.sequence) for record in self.records}
        return len(lengths) == 1

    @property
    def padded_cells(self) -> int:
        """Number of cells added to reach a rectangular matrix (0 for a true alignment)."""

        width = self.length
        return sum(width - len(record.sequence) for record in self.records)

    def padded_sequences(self, pad: str = "?") -> list[str]:
        """Return equal-width sequences.

        Unequal-length input is padded with ``?`` (missing data) rather than ``-`` so that
        artificial padding is never presented as a biological alignment gap.
        """

        width = self.length
        return [record.sequence.ljust(width, pad) for record in self.records]

    def columns(self) -> list[tuple[str, ...]]:
        seqs = self.padded_sequences()
        return [tuple(seq[i] for seq in seqs) for i in range(self.length)]
