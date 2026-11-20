"""Output writer for block, read, mapping, graph and profile-model inputs.

Writes ``summary.json`` (with provenance), ``tables/*.csv``, ``metric_dictionary.json``,
``report.md`` and a small offline ``report.html``. Sequence alignments use the full
:class:`msatk.core.profiler.MSATK` pipeline instead.
"""

from __future__ import annotations

import html
import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from msatk.io import write_csv, write_json
from msatk.report.formatting import format_value
from msatk.stats.registry import metric_dictionary
from msatk.version import __version__

TITLES = {
    "block": "Block (whole-genome) alignment profile",
    "read": "Read alignment QC",
    "mapping": "Mapping / hit table profile",
    "graph": "Sequence graph profile",
    "profile_model": "Profile HMM summary",
}


def write_kind_outputs(
    obj: Any, source: str | Path, outdir: str | Path | None = None, force: bool = False,
    command: str | None = None,
) -> dict[str, Any]:  # fmt: skip
    from msatk.encoding.export import file_sha256

    src = Path(source)
    out = Path(outdir) if outdir else src.with_name(src.name.split(".")[0] + "_msatk")
    if out.exists() and any(out.iterdir()) and not force:
        raise FileExistsError(f"Output directory {out} exists and is not empty; use --force.")
    if out.exists() and force:
        shutil.rmtree(out)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    summary = {
        **obj.summary(),
        "msatk_version": __version__,
        "input": str(src),
        "input_sha256": file_sha256(src) if src.is_file() else None,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "command": command,
        "coordinate_convention": "0-based, half-open (block/read/mapping intervals)",
    }
    tables = obj.tables()
    for name, rows in tables.items():
        write_csv(out / "tables" / f"{name}.csv", rows)
    write_json(out / "summary.json", summary)
    write_json(out / "metric_dictionary.json", metric_dictionary())
    (out / "report.md").write_text(_markdown(summary, tables), encoding="utf-8")
    (out / "report.html").write_text(_html(summary, tables), encoding="utf-8")
    (out / "run_manifest.json").write_text(
        json.dumps(
            {
                "msatk_version": __version__,
                "kind": summary.get("alignment_kind"),
                "files": sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return {"outdir": str(out), "summary": summary, "tables": list(tables)}


_SKIP = {"python_version", "platform", "generated_at", "command", "input_sha256"}


def _markdown(summary: dict[str, Any], tables: dict[str, list[dict[str, Any]]]) -> str:
    kind = str(summary.get("alignment_kind"))
    lines = [f"# {TITLES.get(kind, 'MSATK profile')}", "", f"Input: `{summary.get('input')}`", ""]
    lines += ["| Metric | Value |", "|---|---|"]
    lines += [f"| {k} | {format_value(v, k)} |" for k, v in summary.items() if k not in _SKIP]
    for name, rows in tables.items():
        lines += ["", f"## {name.replace('_', ' ').title()}", ""]
        if not rows:
            lines.append("_No rows._")
            continue
        keys = list(rows[0])
        lines += ["| " + " | ".join(keys) + " |", "|" + "---|" * len(keys)]
        for r in rows[:200]:
            lines.append("| " + " | ".join(format_value(r.get(k), k) for k in keys) + " |")
        if len(rows) > 200:
            lines.append(f"\n_{len(rows) - 200} more rows in tables/{name}.csv._")
    return "\n".join(lines) + "\n"


def _html(summary: dict[str, Any], tables: dict[str, list[dict[str, Any]]]) -> str:
    kind = str(summary.get("alignment_kind"))
    esc = html.escape

    def table(rows: list[dict[str, Any]]) -> str:
        if not rows:
            return "<p class=muted>No rows.</p>"
        keys = list(rows[0])
        head = "".join(f"<th>{esc(k)}</th>" for k in keys)
        body = "".join(
            "<tr>" + "".join(f"<td>{esc(format_value(r.get(k), k))}</td>" for k in keys) + "</tr>"
            for r in rows[:500]
        )
        more = (
            f"<p class=muted>{len(rows) - 500} more rows in the CSV.</p>" if len(rows) > 500 else ""
        )
        return f"<div class=scroll><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>{more}"

    metrics = [{"metric": k, "value": v} for k, v in summary.items() if k not in _SKIP]
    sections = "".join(
        f"<h2>{esc(name.replace('_', ' ').title())}</h2>{table(rows)}"
        for name, rows in tables.items()
    )
    return f"""<!doctype html><html lang=en><head><meta charset=utf-8>
<meta name=viewport content="width=device-width,initial-scale=1">
<title>MSATK {esc(kind)} report</title>
<style>
:root{{--bg:#fff;--fg:#1f2933;--muted:#6b7785;--line:#d9dee3;--head:#f4f6f8}}
@media (prefers-color-scheme:dark){{:root{{--bg:#111418;--fg:#e6e9ec;--muted:#9aa5b1;--line:#2c3238;--head:#1a1f24}}}}
body{{margin:0;padding:24px 16px;background:var(--bg);color:var(--fg);font:14px/1.5 system-ui,sans-serif}}
main{{max-width:1100px;margin:auto}} h1{{font-size:22px}} h2{{font-size:17px;margin-top:28px}}
.muted{{color:var(--muted)}} .scroll{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
th,td{{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left;white-space:nowrap}}
th{{background:var(--head);position:sticky;top:0}}
</style></head><body><main>
<h1>{esc(TITLES.get(kind, "MSATK profile"))}</h1>
<p class=muted>{esc(str(summary.get("input")))} · MSATK {esc(__version__)} · {esc(str(summary.get("generated_at")))}</p>
<h2>Summary</h2>{table(metrics)}{sections}
</main></body></html>
"""
