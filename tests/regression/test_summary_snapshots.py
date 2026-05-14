import json

import pytest

from msatk import profile_alignment


@pytest.mark.regression
def test_tiny_dna_summary_regression(data_dir, tmp_path):
    outdir = tmp_path / "results"
    profile_alignment(
        data_dir / "dna" / "tiny_dna.fasta", outdir=outdir, plots=False, report=False, force=True
    )
    observed = json.loads((outdir / "summary.json").read_text())
    expected = json.loads((data_dir / "dna" / "expected_summary.json").read_text())
    for key, value in expected.items():
        assert observed[key] == value
    assert observed["gap_fraction"] == pytest.approx(1 / 30)
