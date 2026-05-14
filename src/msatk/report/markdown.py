"""Markdown reports."""

from __future__ import annotations


def render_markdown_report(
    summary: dict[str, object],
    warnings: list[str],
    figures: list[str] | None = None,
    interpretation: str = "",
    next_steps: list[str] | None = None,
) -> str:
    figures = figures or []
    next_steps = next_steps or []
    lines = [
        "# MSATK Alignment Report",
        "",
        "One command to profile, visualize, and report multiple sequence alignments.",
        "",
        "## Executive Summary",
        "",
        interpretation,
        "",
        "| Metric | Value |",
        "| --- | --- |",
    ]
    lines.extend(
        f"| {key.replace('_', ' ').title()} | {_format_value(value)} |"
        for key, value in summary.items()
        if key not in {"interpretation", "recommended_next_steps"}
    )
    lines.extend(["", "## QC Warnings", ""])
    lines.extend(f"- {warning}" for warning in warnings)
    lines.extend(["", "## Recommended Next Steps", ""])
    lines.extend(f"- {step}" for step in next_steps)
    lines.extend(["", "## Figures", ""])
    if figures:
        lines.extend(f"![{figure}](figures/{figure})" for figure in figures)
    else:
        lines.append("No figures were generated in this run.")
    lines.append("")
    return "\n".join(lines)


def _format_value(value: object) -> str:
    if isinstance(value, float):
        if 0 <= value <= 1:
            return f"{value:.1%}"
        return f"{value:.6g}"
    return str(value)
