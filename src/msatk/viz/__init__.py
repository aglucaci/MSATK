"""Human-facing visualisation: residue-level alignment figures, logos and themes."""

from __future__ import annotations

from typing import Any

from .themes import THEMES, load_theme

__all__ = ["THEMES", "load_theme", "logo_matrix", "plot_alignment", "plot_logo"]


def __getattr__(name: str) -> Any:  # lazy: avoids importing matplotlib/encoders eagerly
    if name == "plot_alignment":
        from .alignment import plot_alignment

        return plot_alignment
    if name in {"plot_logo", "logo_matrix"}:
        from . import logos

        return getattr(logos, name)
    raise AttributeError(name)
