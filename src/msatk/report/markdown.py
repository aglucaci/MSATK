"""Markdown reports."""

from __future__ import annotations

from msatk.report.formatting import format_value


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
        "![MSATK logo](assets/msatk_logo.png)",
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
    nested = {"specialized", "characterization"}
    lines.extend(
        f"| {key.replace('_', ' ').title()} | {format_value(value, key)} |"
        for key, value in summary.items()
        if key not in {"interpretation", "recommended_next_steps", *nested}
    )
    for key in ("specialized", "characterization"):
        groups = summary.get(key) or {}
        if not isinstance(groups, dict) or not groups:
            continue
        lines.extend(["", f"## {key.title()} analyses", ""])
        for name, values in groups.items():
            lines.extend(
                [
                    f"### {name.replace('_', ' ').capitalize()}",
                    "",
                    "| Metric | Value |",
                    "| --- | --- |",
                ]
            )
            lines.extend(f"| `{k}` | {format_value(v, k)} |" for k, v in (values or {}).items())
            lines.append("")
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
