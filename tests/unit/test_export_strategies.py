"""Export strategies (native / tiles / pad), row ordering and batching."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from msatk.encoding import export_tensor, row_order
from msatk.ml import AlignmentTensorDataset, pad_batch
from msatk.models import Alignment, SequenceRecord


def aln(seqs, ids=None):
    ids = ids or [f"s{i}" for i in range(len(seqs))]
    return Alignment(tuple(SequenceRecord(i, s) for i, s in zip(ids, seqs)))


def manifest(info):
    return json.loads(Path(info["manifest"]).read_text())


def random_dna(n, length, seed=0):
    rng = np.random.default_rng(seed)
    return ["".join(rng.choice(list("ACGT-"), size=length)) for _ in range(n)]


def test_native_keeps_full_shape(tmp_path):
    a = aln(random_dna(5, 17))
    info = export_tensor(a, "dna", tmp_path, encoding="categorical")
    assert info["shape"] == [5, 17] and info["strategy"] == "native"
    assert manifest(info)["strategy"]["alignment_shape"] == [5, 17]
    with pytest.raises(ValueError, match="tiles' and 'pad"):
        export_tensor(a, "dna", tmp_path, strategy="native", width=8)


def test_tiles_cover_every_cell_and_reassemble(tmp_path):
    a = aln(random_dna(7, 23, seed=1))
    full = np.load(export_tensor(a, "dna", tmp_path / "n", encoding="categorical")["tensor"])["x"]
    info = export_tensor(
        a, "dna", tmp_path / "t", encoding="categorical", strategy="tiles", height=3, width=10
    )
    data = np.load(info["tensor"])
    assert data["x"].shape == (3 * 3, 3, 10)
    rebuilt = np.full_like(full, 255)
    for tile, mask, (r0, r1, s0, s1) in zip(data["x"], data["mask"], data["tiles"]):
        h, w = r1 - r0 + 1, s1 - s0 + 1
        assert mask[:h, :w].all() and mask.sum() == h * w
        rebuilt[r0 : r1 + 1, s0 - 1 : s1] = tile[:h, :w]
    assert (rebuilt == full).all()
    m = manifest(info)
    assert m["strategy"]["padded_cells"] == int((~data["mask"]).sum())
    assert m["lossless"] is True


def test_legacy_window_arguments_mean_tiles(tmp_path):
    a = aln(random_dna(4, 20))
    info = export_tensor(a, "dna", tmp_path, window=8, rows_per_tile=2)
    assert info["strategy"] == "tiles"
    assert info["shape"][1:3] == [2, 8]


def test_pad_pads_smaller_and_rejects_larger(tmp_path):
    a = aln(random_dna(3, 10))
    info = export_tensor(
        a, "dna", tmp_path, encoding="categorical", strategy="pad", height=8, width=16
    )
    data = np.load(info["tensor"])
    assert data["x"].shape == (1, 8, 16)
    assert data["mask"][0].sum() == 30
    assert data["tiles"].tolist() == [[0, 2, 1, 10]]
    with pytest.raises(ValueError, match="larger than the pad target"):
        export_tensor(a, "dna", tmp_path, strategy="pad", height=2, width=16)
    with pytest.raises(ValueError, match="needs both"):
        export_tensor(a, "dna", tmp_path, strategy="pad", width=16)


def test_pad_crop_is_recorded_as_lossy(tmp_path):
    a = aln(random_dna(5, 12))
    info = export_tensor(
        a,
        "dna",
        tmp_path,
        encoding="categorical",
        strategy="pad",
        height=4,
        width=10,
        overflow="crop",
    )
    m = manifest(info)
    assert info["shape"] == [1, 4, 10] and m["lossless"] is False
    crop = [t for t in m["transforms"] if t["op"] == "crop_to_pad_target"][0]
    assert crop["rows_dropped"] == 1 and crop["sites_dropped"] == 2


def test_codon_tiles_keep_reading_frame(tmp_path):
    a = aln(["ATGAAACCCGGGTTT", "ATGAAACCAGGGTTA"])
    info = export_tensor(a, "codon", tmp_path, strategy="tiles", width=6, stride=3)
    m = manifest(info)
    assert m["strategy"]["reading_frame"]["preserved"] is True
    assert all((t["site_start"] - 1) % 3 == 0 for t in m["windowing"]["tiles"])
    for bad in ({"width": 7}, {"width": 6, "stride": 4}, {"width": 6, "start": 2}):
        with pytest.raises(ValueError, match="reading frame"):
            export_tensor(a, "codon", tmp_path, strategy="tiles", **bad)
    # DNA alignments have no frame constraint.
    export_tensor(a, "dna", tmp_path, strategy="tiles", width=7)


def test_row_orders():
    a = aln(
        ["AAAAAAAA", "CCCCCCCC", "AAAAAAAC", "CCCCCCCA", "AAAA----"],
        ids=["a1", "c1", "a2", "c2", "a3"],
    )
    ids, d = row_order(a, "dna", "reference", reference="c1")
    assert ids[:2] == ["c1", "c2"] and d["reference"] == "c1"
    ids, d = row_order(a, "dna", "cluster")
    assert sorted(ids) == sorted(a.ids) and d["linkage"] == "average (UPGMA)"
    assert row_order(a, "dna", "cluster")[0] == ids  # deterministic
    ids, d = row_order(a, "dna", "gaps")
    assert ids[-1] == "a3"
    ids, d = row_order(a, "dna", ["c2", "a3"])
    assert ids == ["c2", "a3", "a1", "c1", "a2"] and d["appended_in_input_order"] == 3
    with pytest.raises(ValueError, match="not in the alignment"):
        row_order(a, "dna", ["zz"])
    with pytest.raises(ValueError, match="Unknown row order"):
        row_order(a, "dna", "random")


def test_cluster_groups_similar_sequences_together():
    a = aln(
        ["AAAAAAAAAA", "CCCCCCCCCC", "AAAAAAAAAC", "CCCCCCCCCA", "AAAAAAAACC"],
        ids=["a1", "c1", "a2", "c2", "a3"],
    )
    ids, _ = row_order(a, "dna", "cluster")
    letters = "".join(i[0] for i in ids)
    assert letters in ("aaacc", "ccaaa")


def test_order_file_and_manifest(tmp_path):
    a = aln(random_dna(4, 9), ids=["w", "x", "y", "z"])
    order_file = tmp_path / "order.txt"
    order_file.write_text("# wanted\nz\ny\n")
    info = export_tensor(a, "dna", tmp_path, order=f"file:{order_file}")
    m = manifest(info)
    assert m["sequence_ids"] == ["z", "y", "w", "x"]
    assert m["coordinate_map"]["original_rows_0based"] == [3, 2, 0, 1]
    assert m["row_order"]["method"] == "user"


def test_rgb_encoding_tiles(tmp_path):
    a = aln(["ACGT-", "ACGA-"])
    info = export_tensor(a, "dna", tmp_path, encoding="rgb", strategy="pad", height=4, width=8)
    data = np.load(info["tensor"])
    assert data["x"].shape == (1, 4, 8, 3) and data["x"].dtype == np.uint8
    m = manifest(info)
    assert m["channels"] == ["R", "G", "B"] and m["lossless"] is False


def test_dataset_combines_same_shape_exports_and_rejects_mismatch(tmp_path):
    small = export_tensor(
        aln(random_dna(3, 10)), "dna", tmp_path / "a", strategy="pad", height=4, width=12
    )
    large = export_tensor(
        aln(random_dna(6, 30, seed=2)), "dna", tmp_path / "b", strategy="tiles", height=4, width=12
    )
    ds = AlignmentTensorDataset([small["tensor"], large["tensor"]], channels_first=True)
    assert len(ds) == 1 + 2 * 3
    assert ds[0]["source"] == 0 and ds[1]["source"] == 1
    assert ds[3]["x"].shape == (6, 4, 12)
    other = export_tensor(
        aln(random_dna(3, 10)), "dna", tmp_path / "c", strategy="pad", height=4, width=16
    )
    with pytest.raises(ValueError, match="same --height and --width"):
        AlignmentTensorDataset([small["tensor"], other["tensor"]])
    native = export_tensor(aln(random_dna(3, 10)), "dna", tmp_path / "d")
    with pytest.raises(ValueError, match="Cannot mix"):
        AlignmentTensorDataset([small["tensor"], native["tensor"]])


@pytest.mark.parametrize("channels_first", [False, True])
def test_pad_batch_variable_sizes(tmp_path, channels_first):
    a = export_tensor(aln(random_dna(3, 10)), "dna", tmp_path / "a")
    b = export_tensor(aln(random_dna(5, 7, seed=3)), "dna", tmp_path / "b")
    ds = AlignmentTensorDataset([a["tensor"], b["tensor"]], channels_first=channels_first)
    batch = pad_batch([ds[0], ds[1]], multiple_of=4)
    assert batch["mask"].shape == (2, 8, 12)
    assert batch["mask"][0].sum() == 30 and batch["mask"][1].sum() == 35
    expected = (2, 6, 8, 12) if channels_first else (2, 8, 12, 6)
    assert batch["x"].shape == expected
    assert batch["source"].tolist() == [0, 1]
