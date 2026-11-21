"""End-to-end checks of the v0.2 profile outputs, report and new CLI commands."""

from __future__ import annotations

import json
import re

import pytest

from msatk import MSATK
from msatk.cli import main

mpl = pytest.importorskip("matplotlib")


def test_profile_outputs_and_provenance(data_dir, tmp_path):
    out = tmp_path / "res"
    code = main(
        ["profile", str(data_dir / "codon" / "valid_codon_alignment.fasta"), "--out", str(out),
         "--plots", "png,svg,pdf", "--encodings", "categorical,onehot", "--force"]
    )  # fmt: skip
    assert code == 0
    summary = json.loads((out / "summary.json").read_text())
    assert len(summary["input_sha256"]) == 64
    assert summary["output_schema_version"] == "2.0"
    assert summary["coordinate_convention"] == "1-based, inclusive"
    for name in (
        "metric_dictionary.json",
        "run_manifest.json",
        "qc_findings.json",
        "parameters.yaml",
    ):
        assert (out / name).exists()
    for fmt in ("png", "svg", "pdf"):
        assert (out / "figures" / f"entropy_by_site.{fmt}").exists()
        assert (out / "figures" / f"alignment_overview.{fmt}").exists()
    assert (out / "figures" / "sequence_logo.png").exists()
    assert (out / "tables" / "pairwise_comparisons.csv").exists()
    assert (out / "tables" / "gap_blocks.csv").exists()
    assert (out / "ml" / "alignment_onehot.npz").exists()
    assert not (out / "tables" / "sequence_embeddings.csv").exists()  # embeddings are opt-in
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["config"]["encodings"] == "categorical,onehot"


def test_invalid_plot_format_is_rejected(data_dir, tmp_path):
    assert (
        main(
            [
                "profile",
                str(data_dir / "dna" / "tiny_dna.fasta"),
                "--out",
                str(tmp_path),
                "--plots",
                "gif",
                "--force",
            ]
        )
        == 2
    )


def test_basic_preset_skips_pairwise(data_dir, tmp_path):
    assert (
        main(
            [
                "profile",
                str(data_dir / "dna" / "tiny_dna.fasta"),
                "--out",
                str(tmp_path),
                "--preset",
                "basic",
                "--force",
            ]
        )
        == 0
    )
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["pairwise_method"] == "none"
    assert summary["mean_pairwise_identity"] is None
    assert not (tmp_path / "tables" / "pairwise_identity_matrix.csv").exists()


def test_self_contained_report_has_no_external_references(data_dir, tmp_path):
    profiler = MSATK(data_dir / "protein" / "protein_with_gaps.faa")
    profiler.write_outputs(tmp_path, self_contained=True, force=True)
    html = (tmp_path / "report.html").read_text()
    assert 'src="figures/' not in html and 'href="tables/' not in html
    assert "data:image/png;base64," in html
    assert not re.search(r"""(src|href)=["']https?://""", html)
    data = json.loads(
        re.search(
            r'<script id="msatk-data" type="application/json">(.*?)</script>', html, re.S
        ).group(1)
    )
    assert data["viewer"]["ids"] == profiler.alignment.ids
    assert {"summary", "findings", "sites", "sequences", "metrics"} <= set(data)


def test_report_api_single_file(data_dir, tmp_path):
    path = MSATK(data_dir / "dna" / "tiny_dna.fasta").report(tmp_path / "r.html")
    assert path.exists() and "Alignment viewer" in path.read_text()


@pytest.mark.parametrize(
    "level", ["summary", "site", "sequence", "pairwise", "composition", "gaps"]
)
def test_stats_command(data_dir, tmp_path, level):
    target = tmp_path / f"{level}.json"
    assert (
        main(
            [
                "stats",
                str(data_dir / "dna" / "dna_with_gaps.fasta"),
                "--level",
                level,
                "--output-format",
                "json",
                "--out",
                str(target),
            ]
        )
        == 0
    )
    assert json.loads(target.read_text())


def test_plot_logo_encode_image_commands(data_dir, tmp_path):
    src = str(data_dir / "protein" / "protein_with_gaps.faa")
    assert (
        main(
            [
                "plot",
                src,
                "--out",
                str(tmp_path / "a.svg"),
                "--start",
                "1",
                "--end",
                "5",
                "--tracks",
                "entropy,gaps",
            ]
        )
        == 0
    )
    assert main(["logo", src, "--out", str(tmp_path / "logo.pdf")]) == 0
    assert (
        main(
            [
                "encode",
                src,
                "--out",
                str(tmp_path / "ml"),
                "--encoding",
                "physicochemical",
                "--window",
                "4",
                "--label",
                "family=demo",
            ]
        )
        == 0
    )
    assert main(["image", src, "--out", str(tmp_path / "img.png"), "--scale", "3"]) == 0
    assert main(["image", src, "--out", str(tmp_path / "img224.png"), "--size", "224x224"]) == 0
    for name in ("a.svg", "logo.pdf", "img.png", "img224.png", "img224.png.manifest.json"):
        assert (tmp_path / name).exists()
    assert list((tmp_path / "ml").glob("*.npz"))
    assert main(["metrics", "--format", "json"]) == 0


def test_embeddings_report_actual_method(data_dir):
    rows = MSATK(data_dir / "dna" / "tiny_dna.fasta").embeddings("pca")
    rows = rows.to_dict(orient="records") if hasattr(rows, "to_dict") else rows
    assert {r["method"] for r in rows} <= {"pca", "none (fewer than 3 sequences)"}


def test_report_rebuild(data_dir, tmp_path):
    MSATK(data_dir / "dna" / "tiny_dna.fasta").write_outputs(tmp_path, force=True)
    assert main(["report", str(tmp_path)]) == 0
    assert main(["report", str(tmp_path), "--format", "md"]) == 0
