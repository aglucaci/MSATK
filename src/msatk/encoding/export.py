"""Write ML-ready tensors/images with a self-describing ``manifest.json``."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
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


STRATEGIES = ("native", "tiles", "pad")
OVERFLOW_POLICIES = ("error", "crop")


def export_tensor(
    alignment: Alignment,
    molecule_type: str,
    out_dir: str | Path,
    encoding: str = "onehot",
    ambiguity: str = "missing",
    start: int | None = None,
    end: int | None = None,
    strategy: str | None = None,
    height: int | None = None,
    width: int | None = None,
    stride: int | None = None,
    row_stride: int | None = None,
    overflow: str = "error",
    order: str | Sequence[str] = "input",
    order_reference: str | None = None,
    theme: str | None = "auto",
    window: int | None = None,
    rows_per_tile: int | None = None,
    labels: dict[str, Any] | None = None,
    compress: bool = True,
    stem: str = "alignment",
) -> dict[str, Any]:
    """Encode an alignment and write ``<stem>_<encoding>.npz`` plus a manifest.

    Strategies (no strategy ever resamples the alignment):

    ``native``  the whole alignment at its own size, ``x`` is ``(N, L, C)``.
    ``tiles``   fixed ``height x width`` tiles (default: all rows x ``width``), the
                last tile on each axis padded; ``stride``/``row_stride`` set overlap.
    ``pad``     a single ``height x width`` tensor; smaller alignments are padded.
                Larger ones are rejected unless ``overflow="crop"`` (lossy, recorded).

    ``tiles`` and ``pad`` write ``x`` ``(T, height, width, C)``, ``mask`` (True =
    real cell) and ``tiles`` (``[row_start, row_end, site_start, site_end]`` per
    tile, rows 0-based in the exported row order, sites 1-based in the original
    alignment). ``window``/``rows_per_tile`` are the older names for
    ``width``/``height`` and imply ``strategy="tiles"``.
    """

    if window is not None:
        width = width or window
        strategy = strategy or "tiles"
    if rows_per_tile is not None:
        height = height or rows_per_tile
    strategy = strategy or ("tiles" if width else "native")
    if strategy not in STRATEGIES:
        raise ValueError(f"strategy must be one of {STRATEGIES}, not {strategy!r}.")
    if overflow not in OVERFLOW_POLICIES:
        raise ValueError(f"overflow must be one of {OVERFLOW_POLICIES}, not {overflow!r}.")
    for name, value in (
        ("height", height),
        ("width", width),
        ("stride", stride),
        ("row_stride", row_stride),
    ):
        if value is not None and value < 1:
            raise ValueError(f"{name} must be a positive integer.")
    if strategy == "tiles" and not width:
        raise ValueError("strategy 'tiles' needs a tile width (--width).")
    if strategy == "pad" and not (width and height):
        raise ValueError("strategy 'pad' needs both --height and --width.")
    if strategy == "native" and (width or height):
        raise ValueError("--height/--width apply to the 'tiles' and 'pad' strategies only.")
    frame = _check_reading_frame(molecule_type, strategy, start, width, stride)

    from msatk.encoding.ordering import row_order

    ordered_ids, order_details = row_order(alignment, molecule_type, order, order_reference)
    order_details["computed_on"] = "full alignment (before any site cropping)"
    sub, coords = subset_alignment(
        alignment, start=start, end=end, order=None if order == "input" else ordered_ids
    )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    kwargs: dict[str, Any] = {}
    if encoding == "onehot":
        kwargs["ambiguity"] = ambiguity
    elif encoding == "rgb":
        kwargs["theme"] = theme
    tensor: EncodedTensor = encode(sub, molecule_type, encoding, **kwargs)
    arrays: dict[str, np.ndarray] = {}
    manifest = _base_manifest(sub, coords, molecule_type)
    transforms = list(coords.operations)
    lossless = tensor.lossless
    manifest.update(
        {
            "encoding": encoding,
            "channels": tensor.channels,
            "dtype": str(tensor.array.dtype),
            "encoding_details": tensor.details,
            "row_order": order_details,
            "labels": labels or {},
        }
    )
    n, length = tensor.array.shape[:2]
    strategy_info: dict[str, Any] = {"name": strategy, "alignment_shape": [n, length]}
    if frame:
        strategy_info["reading_frame"] = frame
    if strategy == "native":
        arrays["x"] = tensor.array
        manifest["axes"] = ["sequence", "site", "channel"][: tensor.array.ndim]
    else:
        data = tensor.array
        site_offset = coords.original_sites[0] - 1 if coords.original_sites else 0
        if strategy == "pad":
            assert height is not None and width is not None
            if n > height or length > width:
                if overflow == "error":
                    raise ValueError(
                        f"Alignment is {n} x {length} (sequences x sites), larger than the pad "
                        f"target {height} x {width}. Use strategy 'tiles' to keep every cell, or "
                        "overflow='crop' (--overflow crop) to keep only the first "
                        f"{height} rows and {width} sites."
                    )
                data = data[:height, :width]
                lossless = False
                transforms.append(
                    {
                        "op": "crop_to_pad_target",
                        "rows_dropped": max(n - height, 0),
                        "sites_dropped": max(length - width, 0),
                        "lossless": False,
                    }
                )
            win = sliding_windows(data, width=width, rows=height, site_offset=site_offset)
        else:
            assert width is not None
            win = sliding_windows(
                data,
                width=width,
                stride=stride,
                rows=height,
                row_stride=row_stride,
                site_offset=site_offset,
            )
        arrays["x"] = win.data
        arrays["mask"] = win.mask
        arrays["tiles"] = np.array(
            [[t["row_start"], t["row_end"], t["site_start"], t["site_end"]] for t in win.table],
            dtype=np.int64,
        )
        manifest["axes"] = ["tile", "sequence", "site", "channel"][: win.data.ndim]
        tile_rows, tile_width = int(win.data.shape[1]), int(win.data.shape[2])
        strategy_info.update(
            {
                "height": tile_rows,
                "width": tile_width,
                "stride": (stride or tile_width) if strategy == "tiles" else None,
                "row_stride": (row_stride or tile_rows) if strategy == "tiles" else None,
                "overflow": overflow if strategy == "pad" else None,
                "padded_cells": int((~win.mask).sum()),
            }
        )
        manifest["windowing"] = {
            "width": tile_width,
            "stride": strategy_info["stride"] or tile_width,
            "rows_per_tile": tile_rows,
            "n_tiles": len(win.table),
            "pad_value": 0,
            "mask_semantics": "True = real cell, False = padding",
            "tiles": win.table,
        }
    manifest["strategy"] = strategy_info
    manifest["lossless"] = lossless
    manifest["transforms"] = transforms
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
    return {
        "tensor": str(npz),
        "manifest": str(manifest_path),
        "shape": manifest["shape"],
        "strategy": strategy,
        "lossless": lossless,
    }


def _check_reading_frame(
    molecule_type: str, strategy: str, start: int | None, width: int | None, stride: int | None
) -> dict[str, Any] | None:
    """Codon alignments: tiles must start on codon position 1 and span whole codons."""

    if molecule_type != "codon" or strategy == "native":
        return None
    problems = []
    if start is not None and (start - 1) % 3:
        problems.append(f"start={start} is not codon position 1 (use 1, 4, 7, ...)")
    if width is not None and width % 3:
        problems.append(f"width={width} is not a multiple of 3")
    if stride is not None and stride % 3:
        problems.append(f"stride={stride} is not a multiple of 3")
    if problems:
        raise ValueError(
            "Codon alignment tiles must keep the reading frame: " + "; ".join(problems) + "."
        )
    return {"preserved": True, "tile_first_site_codon_position": 1}


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
