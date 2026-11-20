"""Notebook-friendly public API helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from msatk.core.profiler import MSATK


def profile_alignment(
    alignment: str | Path,
    outdir: str | Path | None = None,
    sequence_type: str = "auto",
    fmt: str = "auto",
    force: bool = False,
    plots: bool = True,
    report: bool = True,
    tables: bool = True,
    validation_mode: str = "permissive",
    pairwise: str = "auto",
    preset: str | None = None,
    plot_formats: str = "png",
    self_contained: bool = False,
    embeddings: str | None = None,
    encodings: list[str] | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Profile an alignment with the same defaults as ``msatk profile``."""

    profiler = MSATK(
        alignment,
        sequence_type=sequence_type,
        fmt=fmt,
        validation_mode=validation_mode,
        pairwise=pairwise,
    )
    return profiler.write_outputs(
        outdir=outdir,
        force=force,
        plots=plots,
        report=report,
        tables=tables,
        preset=preset,
        plot_formats=plot_formats,
        self_contained=self_contained,
        embeddings_method=embeddings,
        encodings=encodings,
        **kwargs,
    )


def as_dataframe(rows: list[dict[str, object]]) -> object:
    """Return a pandas DataFrame when pandas is installed, otherwise return rows unchanged."""

    try:
        import pandas as pd

        return pd.DataFrame(rows)
    except Exception:
        return rows
