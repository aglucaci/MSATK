from msatk.io import read_alignment
from msatk.protein.analysis import amino_acid_composition, residue_class_summary


def test_amino_acid_composition(data_dir):
    alignment = read_alignment(data_dir / "protein" / "protein_with_gaps.faa")
    composition = amino_acid_composition(alignment)
    assert any(row["amino_acid"] == "M" and row["count"] == 3 for row in composition)


def test_residue_class_summary(data_dir):
    alignment = read_alignment(data_dir / "protein" / "protein_with_gaps.faa")
    classes = residue_class_summary(alignment)
    assert any(row["residue_class"] == "hydrophobic" for row in classes)
