"""Write ML-ready tensors/images with a self-describing ``manifest.json``."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from msatk.encoding.images import label_image, resize_nearest, rgb_image, save_image
from msatk.encoding.tensors import EncodedTensor, encode
from msatk.encoding.transforms import CoordinateMap, sliding_windows, subset_alignment
from msatk.models import Alignment
from msatk.version import __version__

MANIFEST_SCHEMA = "msatk.encoding/1.0"


def file_sha256(path: str | Path) -> str | None:
    p = Path(path)
    if not p.is_file():
        return None
    digest = hashlib.sha256()
    with p.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _base_manifest(
    alignment: Alignment, coords: CoordinateMap, molecule_type: str
) -> dict[str, Any]:
    return {
        "schema": MANIFEST_SCHEMA,
        "msatk_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "input": alignment.source,
        "input_sha256": file_sha256(alignment.source) if alignment.source else None,
        "molecule_type": molecule_type,
        "axes": ["sequence", "site", "channel"],
        "coordinate_convention": "sites are 1-based inclusive in the ORIGINAL alignment",
        "sequence_ids": alignment.ids,
        "coordinate_map": coords.as_dict(),
        "padding": {
            "input_padded_cells": alignment.padded_cells,
            "input_padding_symbol": "?",
        },
    }


def export_tensor(
    alignment: Alignment,
    molecule_type: str,
    out_dir: str | Path,
    encoding: str = "onehot",
    ambiguity: str = "missing",
    start: int | None = None,
    end: int | None = None,
    window: int | None = None,
    stride: int | None = None,
    rows_per_tile: int | None = None,
    labels: dict[str, Any] | None = None,
    compress: bool = True,
    stem: str = "alignment",
) -> dict[str, Any]:
    """Encode an alignment (optionally cropped and windowed) and write NPZ + manifest.

    NPZ arrays: ``x`` (tensor), and for windowed exports ``mask`` (real cells) and
    ``tiles`` (``[row_start, row_end, site_start, site_end]`` per tile).
    """

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sub, coords = subset_alignment(alignment, start=start, end=end)
    kwargs = {"ambiguity": ambiguity} if encoding == "onehot" else {}
    tensor: EncodedTensor = encode(sub, molecule_type, encoding, **kwargs)
    arrays: dict[str, np.ndarray] = {}
    manifest = _base_manifest(sub, coords, molecule_type)
    manifest.update(
        {
            "encoding": encoding,
            "channels": tensor.channels,
            "dtype": str(tensor.array.dtype),
            "lossless": tensor.lossless,
            "encoding_details": tensor.details,
            "transforms": coords.operations,
            "labels": labels or {},
        }
    )
    if window:
        site_offset = coords.original_sites[0] - 1 if coords.original_sites else 0
        win = sliding_windows(
            tensor.array, width=window, stride=stride, rows=rows_per_tile, site_offset=site_offset
        )
        arrays["x"] = win.data
        arrays["mask"] = win.mask
        arrays["tiles"] = np.array(
            [[t["row_start"], t["row_end"], t["site_start"], t["site_end"]] for t in win.table],
            dtype=np.int64,
        )
        manifest["axes"] = ["tile", "sequence", "site", "channel"][: win.data.ndim]
        manifest["windowing"] = {
            "width": window,
            "stride": stride or window,
            "rows_per_tile": rows_per_tile or sub.n_sequences,
            "n_tiles": len(win.table),
            "pad_value": 0,
            "mask_semantics": "True = real cell, False = padding",
            "tiles": win.table,
        }
    else:
        arrays["x"] = tensor.array
        if tensor.array.ndim == 2:
            manifest["axes"] = ["sequence", "site"]
    manifest["shape"] = list(arrays["x"].shape)
    npz = out / f"{stem}_{encoding}.npz"
    if compress:
        np.savez_compressed(npz, **arrays)  # type: ignore[arg-type]
    else:
        np.savez(npz, **arrays)  # type: ignore[arg-type]
    manifest["files"] = {"tensor": npz.name}
    manifest_path = out / f"{stem}_{encoding}.manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, default=_json_default), encoding="utf-8"
    )
    return {"tensor": str(npz), "manifest": str(manifest_path), "shape": manifest["shape"]}


def export_image(
    alignment: Alignment,
    molecule_type: str,
    out_path: str | Path,
    kind: str = "rgb",
    theme: str | None = "auto",
    scale: int = 1,
    size: tuple[int, int] | None = None,
    start: int | None = None,
    end: int | None = None,
) -> dict[str, Any]:
    """Write an RGB or label image of the alignment plus a sidecar manifest."""

    sub, coords = subset_alignment(alignment, start=start, end=end)
    target = Path(out_path)
    ops = list(coords.operations)
    if kind == "rgb":
        img = rgb_image(sub, molecule_type, theme=theme, scale=scale)
        legend: Any = {"theme": theme, "note": "RGB colours are arbitrary labels."}
        lossless = size is None
        if scale > 1:
            ops.append({"op": "integer_upscale", "factor": scale, "lossless": True})
    elif kind == "label":
        img, channels = label_image(sub, molecule_type)
        legend = {"pixel_value_to_symbol": dict(enumerate(channels))}
        lossless = size is None
    else:
        raise ValueError("kind must be 'rgb' or 'label'.")
    if size is not None:
        before = img.shape[:2]
        img = resize_nearest(img, size[0], size[1])
        ops.append(
            {
                "op": "resize_nearest",
                "from": list(before),
                "to": list(size),
                "lossless": False,
                "warning": "Resizing merges or duplicates cells; pixel coordinates no longer map "
                "1:1 to alignment cells.",
            }
        )
    save_image(img, target)
    manifest = _base_manifest(sub, coords, molecule_type)
    manifest.update(
        {
            "image_kind": kind,
            "image_shape": list(img.shape),
            "cell_pixels": scale if size is None else None,
            "lossless": lossless,
            "legend": legend,
            "transforms": ops,
            "files": {"image": target.name},
        }
    )
    manifest_path = target.with_suffix(target.suffix + ".manifest.json")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, default=_json_default), encoding="utf-8"
    )
    return {"image": str(target), "manifest": str(manifest_path), "shape": list(img.shape)}


def load_tensor(path: str | Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Load an exported NPZ and its manifest."""

    npz_path = Path(path)
    manifest_path = npz_path.with_name(npz_path.name.replace(".npz", ".manifest.json"))
    with np.load(npz_path) as data:
        arrays = {k: data[k] for k in data.files}
    manifest = (
        json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    )
    return arrays, manifest


def _json_default(value: object) -> object:
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)
