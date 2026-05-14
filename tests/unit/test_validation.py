import pytest

from msatk.exceptions import AlignmentFormatError, AlignmentValidationError
from msatk.io import read_alignment


def test_empty_file_raises_clear_error(data_dir):
    with pytest.raises(AlignmentFormatError, match="empty"):
        read_alignment(data_dir / "edge_cases" / "empty_file.fasta")


def test_duplicate_ids_fail_in_strict_mode(data_dir):
    with pytest.raises(AlignmentValidationError, match="Duplicate"):
        read_alignment(data_dir / "edge_cases" / "duplicate_ids.fasta", validation_mode="strict")


def test_unequal_lengths_fail_in_strict_mode(data_dir):
    with pytest.raises(AlignmentValidationError, match="unequal"):
        read_alignment(data_dir / "dna" / "unequal_lengths.fasta", validation_mode="strict")
