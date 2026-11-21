"""ML encodings: reversibility, channel semantics, windowing and manifests."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from msatk.encoding import (
    categorical,
    decode_categorical,
    export_image,
    export_tensor,
    feature_channels,
    full_vocabulary,
    onehot,
    physicochemical,
    rgb_image,
    sliding_windows,
    subset_alignment,
)
from msatk.ml import AlignmentTensorDataset
from msatk.models import Alignment, SequenceRecord


def aln(*seqs):
    return Alignment(tuple(SequenceRecord(f"s{i}", s) for i, s in enumerate(seqs)))


@pytest.mark.parametrize(
    "seqs,kind",
    [
        (("ACGT-N", "ACRY.?", "acgt-n"), "dna"),
        (("MKV-XB*", "MNVQ.Z?"), "protein"),
    ],
)
def test_categorical_round_trip_is_exact(seqs, kind):
    a = aln(*seqs)
    t = categorical(a, kind)
    assert t.lossless
    assert t.array.dtype == np.uint8 and t.shape == (len(seqs), len(seqs[0]))
    back = decode_categorical(t.array, t.ids, full_vocabulary(kind))
    assert back.sequences == [s.upper() for s in seqs]


def test_categorical_flags_unknown_symbols():
    t = categorical(aln("AC#T"), "dna")
    assert not t.lossless and t.details["other_cells"] == 1


@pytest.mark.parametrize("policy", ["missing", "fractional", "separate"])
def test_onehot_cells_sum_to_one(policy):
    t = onehot(aln("ACGTRN-?", "ACGTYB.-"), "dna", ambiguity=policy)
    assert np.allclose(t.array.sum(axis=2), 1.0)
    assert t.channels[:5] == ["A", "C", "G", "T", "<gap>"]


def test_onehot_fractional_iupac():
    t = onehot(aln("RN?"), "dna", ambiguity="fractional")
    r, n, q = t.array[0]
    assert r[:4].tolist() == [0.5, 0.0, 0.5, 0.0]
    assert n[:4].tolist() == [0.25] * 4
    assert q[-1] == 1.0 and q[:4].sum() == 0


def test_protein_channels():
    p = physicochemical(aln("IK-X"), "protein")
    assert p.shape == (1, 4, 7)
    assert p.array[0, 0, 0] == pytest.approx(1.0)  # I hydropathy 4.5/4.5
    assert p.array[0, 1, 1] == pytest.approx(1.0)  # K charge
    assert p.array[0, 2, 5] == 1.0 and p.array[0, 3, 6] == 1.0
    with pytest.raises(ValueError):
        physicochemical(aln("ACGT"), "dna")
    f = feature_channels(aln("ACGT", "ACGT", "ACGA"), "dna")
    assert f.channels[4] == "matches_consensus"
    assert f.array[2, 3, 4] == 0.0 and f.array[0, 3, 4] == 1.0


def test_sliding_windows_keep_coordinates_and_masks():
    x = np.arange(2 * 10).reshape(2, 10)
    w = sliding_windows(x, width=4, stride=3, site_offset=100)
    assert w.data.shape == (3, 2, 4)
    assert [(t["site_start"], t["site_end"]) for t in w.table] == [
        (101, 104),
        (104, 107),
        (107, 110),
    ]
    w2 = sliding_windows(x, width=4)
    assert w2.data.shape == (3, 2, 4)
    assert w2.mask[-1].sum() == 2 * 2  # last window has 2 real sites
    assert w2.table[-1]["valid_sites"] == 2
    assert (w2.data[-1][~w2.mask[-1]] == 0).all()


def test_subset_alignment_coordinate_map():
    a = aln("AC-GT", "ACGGT", "A--GT")
    sub, cmap = subset_alignment(
        a, start=2, end=5, ids=["s2", "s0"], order=["s2", "s0"], max_site_gap_fraction=0.0
    )
    assert sub.ids == ["s2", "s0"]
    # site 2 is gapped in s2 and site 3 in both rows, so only sites 4-5 survive the gap filter
    assert cmap.original_sites == [4, 5]
    assert sub.sequences == ["GT", "GT"]
    assert [op["op"] for op in cmap.operations] == [
        "crop_sites",
        "select_rows",
        "reorder_rows",
        "drop_gappy_sites",
    ]


def test_export_tensor_and_dataset(tmp_path):
    a = aln("ACGTACGTAC", "ACGTTCGTAC")
    info = export_tensor(
        a, "dna", tmp_path, encoding="onehot", window=4, stride=4, labels={"split": "train"}
    )
    manifest = json.loads(Path(info["manifest"]).read_text())
    assert manifest["schema"].startswith("msatk.encoding/")
    assert manifest["windowing"]["n_tiles"] == 3
    assert manifest["labels"] == {"split": "train"}
    assert manifest["channels"][-1] == "<missing>"
    ds = AlignmentTensorDataset(info["tensor"], channels_first=True)
    assert len(ds) == 3
    item = ds[2]
    assert item["x"].shape == (6, 2, 4)
    assert item["tile"].tolist() == [0, 1, 9, 10]
    assert item["mask"].sum() == 4


def test_export_unwindowed_categorical(tmp_path):
    info = export_tensor(
        aln("ACGT", "AC-T"), "dna", tmp_path, encoding="categorical", start=2, end=3
    )
    manifest = json.loads(Path(info["manifest"]).read_text())
    assert manifest["lossless"] is True
    assert manifest["coordinate_map"]["original_sites_1based"] == [2, 3]
    assert info["shape"] == [2, 2]


def test_images(tmp_path):
    a = aln("ACGT", "A-GT")
    img = rgb_image(a, "dna", scale=3)
    assert img.shape == (6, 12, 3) and img.dtype == np.uint8
    info = export_image(a, "dna", tmp_path / "a.png", scale=2)
    assert (tmp_path / "a.png").exists()
    m = json.loads(Path(info["manifest"]).read_text())
    assert m["lossless"] is True and m["cell_pixels"] == 2
    info = export_image(a, "dna", tmp_path / "b.png", size=(8, 8))
    m = json.loads(Path(info["manifest"]).read_text())
    assert m["lossless"] is False and m["transforms"][-1]["op"] == "resize_nearest"
    info = export_image(a, "dna", tmp_path / "c.png", kind="label")
    m = json.loads(Path(info["manifest"]).read_text())
    assert m["legend"]["pixel_value_to_symbol"]["16"] == "-"


def test_custom_theme_file(tmp_path):
    theme = tmp_path / "theme.yaml"
    theme.write_text('A: "#000000"\n# comment\nC: "#ffffff"\n')
    img = rgb_image(aln("AC"), "dna", theme=str(theme))
    assert img[0, 0].tolist() == [0, 0, 0] and img[0, 1].tolist() == [255, 255, 255]
