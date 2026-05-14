import importlib.util
from pathlib import Path

import pytest

from msatk.io.readers import read_alignment


@pytest.mark.skipif(importlib.util.find_spec("pysam") is None, reason="pysam is required")
def test_read_tiny_bam_fixture(data_dir):
    if not (data_dir / "formats" / "tiny.bam").exists():
        pytest.skip("tiny.bam has not been generated")
    alignment = read_alignment(data_dir / "formats" / "tiny.bam")
    assert alignment.n_sequences == 3


@pytest.mark.skipif(importlib.util.find_spec("pysam") is None, reason="pysam is required")
def test_read_tiny_cram_fixture(data_dir):
    if not (data_dir / "formats" / "tiny.cram").exists():
        pytest.skip("tiny.cram has not been generated")
    alignment = read_alignment(data_dir / "formats" / "tiny.cram")
    assert alignment.n_sequences == 3


def test_binary_fixture_generators_are_present(data_dir):
    formats = Path(data_dir) / "formats"
    assert (formats / "tiny.bam").exists()
    assert (formats / "reference.fasta").exists()
    assert (formats / "create_tiny_bam.py").exists()
    assert (formats / "create_binary_alignment_fixtures.py").exists()
