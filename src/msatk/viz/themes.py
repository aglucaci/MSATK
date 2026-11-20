"""Residue colour themes for human-facing alignment images.

RGB colours are arbitrary *labels* chosen for legibility; they carry no
biochemical distance information and must not be used as ML features. Use
:mod:`msatk.encoding.tensors` for scientifically faithful encodings.
"""

from __future__ import annotations

import json
from pathlib import Path

GAP_COLOR = "#f1f3f5"
MISSING_COLOR = "#ced4da"
OTHER_COLOR = "#868e96"

THEMES: dict[str, dict[str, str]] = {
    "nucleotide": {
        "A": "#3fa34d", "C": "#2f6fdf", "G": "#f2a900", "T": "#d7263d", "U": "#d7263d",
        "R": "#9bbf85", "Y": "#9aa9e0", "S": "#a7c3e8", "W": "#e3a3a3", "K": "#e8c98a",
        "M": "#a8d5a2", "B": "#bfbfbf", "D": "#bfbfbf", "H": "#bfbfbf", "V": "#bfbfbf",
        "N": MISSING_COLOR,
    },
    # Physicochemical grouping: hydrophobic, aromatic, positive, negative, polar, special.
    "protein": {
        "A": "#80a0f0", "I": "#80a0f0", "L": "#80a0f0", "M": "#80a0f0", "V": "#80a0f0",
        "F": "#5a6fd6", "W": "#5a6fd6", "Y": "#15a4a4",
        "K": "#f01505", "R": "#f01505", "H": "#e57c2f",
        "D": "#c048c0", "E": "#c048c0",
        "S": "#15c015", "T": "#15c015", "N": "#15c015", "Q": "#15c015",
        "C": "#f08080", "G": "#f09048", "P": "#c0c000",
        "B": "#bfbfbf", "Z": "#bfbfbf", "J": "#bfbfbf", "X": MISSING_COLOR, "*": "#000000",
    },
    # Colour-vision-deficiency-safe nucleotide palette (Okabe–Ito).
    "nucleotide_cvd": {
        "A": "#009e73", "C": "#0072b2", "G": "#e69f00", "T": "#d55e00", "U": "#d55e00",
        "N": MISSING_COLOR,
    },
    "grayscale": {},
}  # fmt: skip


def hex_to_rgb(color: str) -> tuple[int, int, int]:
    c = color.lstrip("#")
    return int(c[0:2], 16), int(c[2:4], 16), int(c[4:6], 16)


def load_theme(
    theme: str | dict[str, str] | Path | None, molecule_type: str = "dna"
) -> dict[str, str]:
    """Resolve a theme name, mapping, or JSON/YAML file into a symbol→hex mapping.

    Gap (``-``/``.``) and missing (``?``) colours are always defined and may be overridden.
    """

    if theme is None or theme == "auto":
        theme = "protein" if molecule_type in {"protein", "translated_cds"} else "nucleotide"
    if isinstance(theme, dict):
        mapping = dict(theme)
    elif isinstance(theme, str) and theme in THEMES:
        mapping = dict(THEMES[theme])
    else:
        path = Path(theme)
        if not path.exists():
            raise ValueError(
                f"Unknown theme {theme!r}; choose from {sorted(THEMES)} or a JSON/YAML file."
            )
        text = path.read_text(encoding="utf-8")
        if path.suffix.lower() in {".yaml", ".yml"}:
            try:
                import yaml  # type: ignore
            except ImportError:
                mapping = _simple_yaml(text)
            else:
                mapping = {str(k): str(v) for k, v in (yaml.safe_load(text) or {}).items()}
        else:
            mapping = {str(k): str(v) for k, v in json.loads(text).items()}
    mapping = {k.upper() if len(k) == 1 else k: v for k, v in mapping.items()}
    mapping.setdefault("-", GAP_COLOR)
    mapping.setdefault(".", mapping["-"])
    mapping.setdefault("?", MISSING_COLOR)
    return mapping


def _simple_yaml(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        out[key.strip().strip("'\"")] = value.strip().strip("'\"")
    return out
