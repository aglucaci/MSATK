"""Shared value formatting for reports."""

from __future__ import annotations

PERCENT_HINTS = ("fraction", "identity", "occupancy", "gc_content", "support")


def format_value(value: object, key: str = "") -> str:
    """Format a summary value; fractions as percentages only when the field is a proportion."""

    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        if any(hint in key for hint in PERCENT_HINTS) and "skew" not in key and 0 <= value <= 1:
            return f"{value:.1%}"
        return f"{value:.4g}"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, list):
        return "; ".join(str(v) for v in value)
    if isinstance(value, dict):
        return "; ".join(
            f"{k} {v:.4g}" if isinstance(v, float) else f"{k} {format_value(v, str(k))}"
            for k, v in value.items()
        )
    return str(value)


def label(key: str) -> str:
    return key.replace("_", " ").capitalize()
