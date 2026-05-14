"""Standalone HTML reports."""

from __future__ import annotations

from html import escape


def render_html_report(
    summary: dict[str, object],
    warnings: list[str],
    figures: list[str] | None = None,
    interpretation: str = "",
    next_steps: list[str] | None = None,
    title: str = "MSATK Alignment Report",
) -> str:
    figures = figures or []
    next_steps = next_steps or []
    metric_keys = [
        "number_of_sequences",
        "alignment_length",
        "molecule_type",
        "gap_fraction",
        "variable_sites",
        "mean_pairwise_identity",
    ]
    cards = "\n".join(
        f'<div class="metric"><span>{escape(key.replace("_", " ").title())}</span><strong>{escape(_format_value(summary.get(key, "")))}</strong></div>'
        for key in metric_keys
    )
    summary_rows = "\n".join(
        f"<tr><th>{escape(str(key).replace('_', ' ').title())}</th><td>{escape(_format_value(value))}</td></tr>"
        for key, value in summary.items()
        if key not in {"interpretation", "recommended_next_steps"}
    )
    warning_items = "\n".join(f"<li>{escape(item)}</li>" for item in warnings)
    step_items = "\n".join(f"<li>{escape(step)}</li>" for step in next_steps)
    figure_blocks = "\n".join(
        f'<figure><img src="figures/{escape(path)}" alt="{escape(path)}"><figcaption>{escape(path)}</figcaption></figure>'
        for path in figures
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>
    :root {{ color-scheme: light; --ink: #17202a; --muted: #5b6773; --line: #d8dee6; --brand: #0f766e; --wash: #f6faf9; --alert: #9a3412; }}
    body {{ margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: var(--ink); background: #fff; }}
    header {{ padding: 34px 40px; background: var(--wash); border-bottom: 1px solid var(--line); }}
    main {{ max-width: 1160px; margin: 0 auto; padding: 28px 24px 48px; }}
    h1 {{ margin: 0 0 6px; font-size: 34px; letter-spacing: 0; }}
    h2 {{ margin-top: 34px; font-size: 20px; }}
    a {{ color: var(--brand); }}
    .brand {{ color: var(--brand); font-weight: 800; text-transform: uppercase; letter-spacing: 0; }}
    .subtitle, .muted {{ color: var(--muted); }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; margin: 24px 0; }}
    .metric {{ border: 1px solid var(--line); padding: 14px; border-radius: 8px; background: #fff; }}
    .metric span {{ display: block; color: var(--muted); font-size: 13px; }}
    .metric strong {{ display: block; margin-top: 6px; font-size: 22px; }}
    .callout {{ border-left: 4px solid var(--brand); background: var(--wash); padding: 14px 16px; }}
    table {{ width: 100%; border-collapse: collapse; border: 1px solid var(--line); }}
    th, td {{ text-align: left; padding: 9px 11px; border-bottom: 1px solid var(--line); font-size: 14px; vertical-align: top; }}
    th {{ width: 34%; background: #fbfcfd; }}
    ul {{ padding-left: 22px; }}
    .warnings li {{ color: var(--alert); }}
    figure {{ margin: 18px 0; border: 1px solid var(--line); padding: 12px; border-radius: 8px; }}
    img {{ max-width: 100%; height: auto; display: block; }}
    figcaption {{ color: var(--muted); margin-top: 8px; font-size: 13px; }}
    .links {{ display: flex; flex-wrap: wrap; gap: 12px; }}
    .links a {{ border: 1px solid var(--line); border-radius: 8px; padding: 8px 10px; text-decoration: none; }}
    footer {{ color: var(--muted); margin-top: 36px; font-size: 13px; }}
  </style>
</head>
<body>
  <header>
    <div class="brand">MSATK</div>
    <h1>{escape(title)}</h1>
    <p class="subtitle">One command to profile, visualize, and report multiple sequence alignments.</p>
  </header>
  <main>
    <section>
      <h2>Executive Summary</h2>
      <div class="metrics">{cards}</div>
      <p class="callout">{escape(interpretation)}</p>
    </section>
    <section>
      <h2>Input And Detection</h2>
      <table>{summary_rows}</table>
    </section>
    <section>
      <h2>Alignment Quality</h2>
      <p class="muted">Review high-gap sequences, high-gap sites, and any outlier warnings before downstream analyses.</p>
      <div class="links">
        <a href="tables/per_sequence_stats.csv">Per-sequence statistics</a>
        <a href="tables/per_site_stats.csv">Per-site statistics</a>
        <a href="tables/qc_flags.csv">QC flags</a>
      </div>
    </section>
    <section>
      <h2>Composition</h2>
      <div class="links">
        <a href="tables/composition_summary.csv">Composition table</a>
        <a href="tables/alignment_summary.csv">Alignment summary</a>
      </div>
    </section>
    <section>
      <h2>Diversity And Conservation</h2>
      <div class="links">
        <a href="tables/entropy_by_site.csv">Entropy by site</a>
        <a href="tables/pairwise_identity_matrix.csv">Pairwise identity matrix</a>
        <a href="tables/distance_matrix.csv">Distance matrix</a>
      </div>
    </section>
    <section>
      <h2>Figures</h2>
      {figure_blocks if figure_blocks else "<p>No figures were generated in this run. Install the plots extra to enable PNG output.</p>"}
    </section>
    <section>
      <h2>QC Warnings</h2>
      <ul class="warnings">{warning_items}</ul>
    </section>
    <section>
      <h2>Recommended Next Steps</h2>
      <ul>{step_items}</ul>
    </section>
    <section>
      <h2>Methods And Parameters</h2>
      <p>Generated by MSATK. See <a href="parameters.yaml">parameters.yaml</a>, <a href="summary.json">summary.json</a>, and <a href="msatk.log">msatk.log</a>.</p>
    </section>
    <footer>Generated by MSATK.</footer>
  </main>
</body>
</html>
"""


def _format_value(value: object) -> str:
    if isinstance(value, float):
        if 0 <= value <= 1:
            return f"{value:.1%}"
        return f"{value:.6g}"
    return str(value)
