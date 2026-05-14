"""Protein-specific Python API."""

from __future__ import annotations

from pathlib import Path

from msatk.core.profiler import AlignmentProfiler
from msatk.io import write_csv
from msatk.protein.analysis import (
    amino_acid_composition,
    hydrophobicity_summary,
    residue_class_summary,
)


class ProteinProfiler(AlignmentProfiler):
    """Protein alignment profiler."""

    def __init__(self, path: str | Path, sequence_type: str = "protein", fmt: str = "auto") -> None:
        super().__init__(path, sequence_type=sequence_type, fmt=fmt)

    def amino_acid_composition(self) -> list[dict[str, object]]:
        return amino_acid_composition(self.alignment)

    def residue_class_summary(self) -> list[dict[str, object]]:
        return residue_class_summary(self.alignment)

    def hydrophobicity_summary(self) -> list[dict[str, object]]:
        return hydrophobicity_summary(self.alignment)

    def write_protein_tables(self, out_dir: str | Path) -> None:
        out = Path(out_dir)
        write_csv(out / "amino_acid_summary.csv", self.amino_acid_composition())
        write_csv(out / "residue_class_summary.csv", self.residue_class_summary())
        write_csv(out / "hydrophobicity_summary.csv", self.hydrophobicity_summary())
