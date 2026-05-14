"""Codon-aware Python API."""

from __future__ import annotations

from pathlib import Path

from msatk.codon.analysis import (
    amino_acid_from_codons,
    codon_position_entropy,
    codon_position_gc,
    codon_usage,
    rscu,
    stop_codon_report,
)
from msatk.core.profiler import AlignmentProfiler
from msatk.io import write_csv


class CodonProfiler(AlignmentProfiler):
    """Codon-aware alignment profiler."""

    def __init__(self, path: str | Path, sequence_type: str = "codon", fmt: str = "auto") -> None:
        super().__init__(path, sequence_type=sequence_type, fmt=fmt)

    def codon_usage(self) -> list[dict[str, object]]:
        return codon_usage(self.alignment)

    def rscu(self) -> list[dict[str, object]]:
        return rscu(self.alignment)

    def amino_acid_from_codons(self) -> list[dict[str, object]]:
        return amino_acid_from_codons(self.alignment)

    def stop_codon_report(self) -> list[dict[str, object]]:
        return stop_codon_report(self.alignment)

    def codon_position_gc(self) -> list[dict[str, object]]:
        return codon_position_gc(self.alignment)

    def codon_position_entropy(self) -> list[dict[str, object]]:
        return codon_position_entropy(self.per_site_stats())

    def write_codon_tables(self, out_dir: str | Path) -> None:
        out = Path(out_dir)
        write_csv(out / "codon_usage.csv", self.codon_usage())
        write_csv(out / "rscu.csv", self.rscu())
        write_csv(out / "amino_acid_from_codons.csv", self.amino_acid_from_codons())
        write_csv(out / "stop_codon_report.csv", self.stop_codon_report())
        write_csv(out / "codon_position_gc.csv", self.codon_position_gc())
        write_csv(out / "codon_position_entropy.csv", self.codon_position_entropy())
