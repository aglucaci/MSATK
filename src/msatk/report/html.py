"""Offline HTML reports.

``render_interactive_report`` produces a dashboard that works from the local
filesystem with no server, account or network request: all CSS/JS is inline and
all data is embedded as JSON. It contains a scrollable/zoomable residue viewer,
interactive site tracks, sortable/searchable tables with CSV download, a
pairwise heatmap, QC findings and full provenance. With ``self_contained=True``
figures are embedded as data URIs so the report is a single portable file.
"""

from __future__ import annotations

import base64
import json
import math
from html import escape
from pathlib import Path
from typing import Any

from msatk.core.qc import QCFinding
from msatk.report.formatting import format_value, label
from msatk.stats.registry import metric_dictionary
from msatk.version import __version__

DEFAULT_MAX_VIEWER_CELLS = 3_000_000
MAX_HEATMAP_SEQUENCES = 400

SITE_COLUMNS = [
    "site_index", "site_class", "occupancy", "gap_fraction", "dominant_character",
    "conservation_score", "consensus_support", "entropy", "normalized_entropy",
    "observed_state_count", "parsimony_informative", "codon_position", "gc_content",
]  # fmt: skip
SEQUENCE_COLUMNS = [
    "row_index", "sequence_id", "ungapped_length", "span_start", "span_end", "gap_fraction",
    "internal_gap_count", "max_gap_run", "ambiguous_count", "missing_count", "occupancy",
    "consensus_identity", "consensus_identity_robust_z", "composition_distance",
    "composition_robust_z", "mean_identity_to_others", "nearest_neighbor",
    "nearest_neighbor_identity", "gc_content", "duplicate_of",
]  # fmt: skip
CARD_KEYS = [
    ("number_of_sequences", "Sequences"),
    ("alignment_length", "Sites"),
    ("molecule_type", "Type"),
    ("occupancy", "Occupancy"),
    ("variable_sites", "Variable sites"),
    ("parsimony_informative_sites", "Parsimony-informative"),
    ("mean_pairwise_identity", "Mean pairwise identity"),
    ("entropy_mean", "Mean entropy (bits)"),
]


def _clean(value: Any) -> Any:
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return None
        return round(value, 6)
    return value


def _rows(rows: list[dict[str, Any]], columns: list[str]) -> dict[str, Any]:
    cols = [c for c in columns if any(c in r for r in rows[:1])] if rows else []
    return {"columns": cols, "rows": [[_clean(r.get(c)) for c in cols] for r in rows]}


def _json_for_script(data: Any) -> str:
    return (
        json.dumps(data, separators=(",", ":"), default=str)
        .replace("</", "<\\/")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def _figure_src(base_dir: Path | None, name: str, self_contained: bool) -> str | None:
    path = (base_dir / "figures" / name) if base_dir else None
    if not self_contained:
        return f"figures/{name}"
    if path is None or not path.exists():
        return None
    mime = {"png": "image/png", "svg": "image/svg+xml", "jpg": "image/jpeg"}.get(
        path.suffix.lstrip(".").lower()
    )
    if mime is None:
        return None
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def render_interactive_report(
    summary: dict[str, Any],
    bundle: Any,
    findings: list[QCFinding],
    interpretation_text: str = "",
    next_steps: list[str] | None = None,
    figures: list[str] | None = None,
    base_dir: Path | None = None,
    self_contained: bool = False,
    embeddings: list[dict[str, Any]] | None = None,
    max_viewer_cells: int = DEFAULT_MAX_VIEWER_CELLS,
    title: str = "MSATK Alignment Report",
) -> str:
    """Render the interactive offline report as one HTML string."""

    from msatk.viz.themes import load_theme

    enc = bundle.encoded
    alignment = enc.alignment
    cells = alignment.n_sequences * alignment.length
    viewer: dict[str, Any] | None = None
    if cells <= max_viewer_cells:
        viewer = {
            "ids": alignment.ids,
            "rows": [s.upper() for s in alignment.padded_sequences()],
            "consensus": enc.consensus_string(),
        }
    pw = bundle.pairwise
    heatmap = None
    if pw is not None and pw.is_exact and len(pw.ids) <= MAX_HEATMAP_SEQUENCES:
        heatmap = {
            "ids": pw.ids,
            "values": [
                [None if v != v else round(float(v), 4) for v in row] for row in pw.identity
            ],
        }
    figure_list = [
        {"name": f, "src": _figure_src(base_dir, f, self_contained)}
        for f in (figures or [])
        if f.rsplit(".", 1)[-1].lower() in {"png", "svg", "jpg"}
    ]
    figure_list = [f for f in figure_list if f["src"]]
    data = {
        "version": __version__,
        "summary": {k: _clean(v) for k, v in summary.items()},
        "findings": [f.as_dict() | {"items": f.items[:500]} for f in findings],
        "sites": _rows(bundle.sites, SITE_COLUMNS),
        "sequences": _rows(bundle.sequences, SEQUENCE_COLUMNS),
        "viewer": viewer,
        "viewerLimit": max_viewer_cells,
        "theme": load_theme("auto", enc.policy.molecule_type),
        "heatmap": heatmap,
        "pairSummary": None if pw is None else pw.summary(),
        "embeddings": embeddings or [],
        "figures": figure_list,
        "metrics": metric_dictionary(),
        "selfContained": self_contained,
    }
    cards = "\n".join(
        f'<div class="card"><span>{escape(lab)}</span><strong>'
        f"{escape(format_value(summary.get(key), key))}</strong></div>"
        for key, lab in CARD_KEYS
    )
    steps = "".join(f"<li>{escape(s)}</li>" for s in (next_steps or []))
    logo = (
        '<div class="wordmark" aria-hidden="true">MSA<b>TK</b></div>'
        if self_contained
        else '<img class="logo" src="assets/msatk_logo.png" alt="MSATK logo" '
        "onerror=\"this.replaceWith(Object.assign(document.createElement('div'),"
        "{className:'wordmark',innerHTML:'MSA<b>TK</b>'}))\">"
    )
    links = "" if self_contained else _download_links()
    return (
        _TEMPLATE.replace("{{TITLE}}", escape(title))
        .replace("{{LOGO}}", logo)
        .replace("{{INPUT}}", escape(str(summary.get("input", ""))))
        .replace("{{GENERATED}}", escape(str(summary.get("generated_at", ""))))
        .replace("{{VERSION}}", escape(__version__))
        .replace("{{CARDS}}", cards)
        .replace("{{INTERPRETATION}}", escape(interpretation_text))
        .replace("{{STEPS}}", steps)
        .replace("{{LINKS}}", links)
        .replace("{{DATA}}", _json_for_script(data))
    )


def _download_links() -> str:
    files = [
        ("summary.json", "Summary (JSON)"),
        ("tables/per_site_stats.csv", "Per-site statistics"),
        ("tables/per_sequence_stats.csv", "Per-sequence statistics"),
        ("tables/pairwise_comparisons.csv", "Pairwise comparisons"),
        ("tables/gap_blocks.csv", "Gap blocks"),
        ("tables/qc_findings.csv", "QC findings"),
        ("metric_dictionary.json", "Metric dictionary"),
        ("parameters.yaml", "Parameters"),
        ("run_manifest.json", "Run manifest"),
    ]
    items = "".join(f'<a href="{href}">{escape(text)}</a>' for href, text in files)
    return f'<h3>Output files</h3><div class="links">{items}</div>'


def render_html_report(
    summary: dict[str, object],
    warnings: list[str],
    figures: list[str] | None = None,
    interpretation: str = "",
    next_steps: list[str] | None = None,
    title: str = "MSATK Alignment Report",
) -> str:
    """Static summary report rebuilt from an existing output directory (``msatk report``)."""

    rows = "\n".join(
        f"<tr><th>{escape(label(str(k)))}</th><td>{escape(format_value(v, str(k)))}</td></tr>"
        for k, v in summary.items()
        if k not in {"interpretation", "recommended_next_steps"}
    )
    warn = "".join(f"<li>{escape(w)}</li>" for w in warnings)
    steps = "".join(f"<li>{escape(s)}</li>" for s in (next_steps or []))
    figs = "".join(
        f'<figure><img src="figures/{escape(f)}" alt="{escape(f)}"><figcaption>{escape(f)}</figcaption></figure>'
        for f in (figures or [])
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<style>{_BASE_CSS}
th{{width:34%}} figure img{{max-width:100%}}</style></head>
<body><header class="top"><div class="wrap"><div class="wordmark">MSA<b>TK</b></div><h1>{escape(title)}</h1></div></header>
<main class="wrap">
<section><h2>Summary</h2><p class="callout">{escape(interpretation)}</p>
<table class="kv">{rows}</table></section>
<section><h2>QC findings</h2><ul>{warn}</ul></section>
<section><h2>Recommended next steps</h2><ul>{steps}</ul></section>
<section><h2>Figures</h2>{figs or "<p>No figures.</p>"}</section>
<footer>Generated by MSATK {escape(__version__)}.</footer></main></body></html>
"""


_BASE_CSS = """
:root{--bg:#ffffff;--panel:#f7f9fa;--ink:#17202a;--muted:#5b6773;--line:#dde3e8;--brand:#0f766e;--brand-soft:#e6f3f1;--warn:#b45309;--err:#b91c1c;--info:#1d4ed8;--ok:#15803d;color-scheme:light}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#0f1417;--panel:#161d21;--ink:#e6ecef;--muted:#98a6b0;--line:#2a343a;--brand:#2dd4bf;--brand-soft:#123430;--warn:#f59e0b;--err:#f87171;--info:#60a5fa;--ok:#4ade80;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#0f1417;--panel:#161d21;--ink:#e6ecef;--muted:#98a6b0;--line:#2a343a;--brand:#2dd4bf;--brand-soft:#123430;--warn:#f59e0b;--err:#f87171;--info:#60a5fa;--ok:#4ade80;color-scheme:dark}
*{box-sizing:border-box}
body{margin:0;font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;color:var(--ink);background:var(--bg)}
.wrap{max-width:1240px;margin:0 auto;padding:0 16px}
header.top{background:var(--panel);border-bottom:1px solid var(--line);padding:22px 0}
h1{margin:4px 0;font-size:26px} h2{font-size:19px;margin:0 0 12px} h3{font-size:15px;margin:18px 0 8px}
section{padding:26px 0;border-bottom:1px solid var(--line)}
.wordmark{font-weight:800;letter-spacing:.04em;color:var(--brand);font-size:20px}.wordmark b{color:var(--ink)}
.callout{border-left:4px solid var(--brand);background:var(--brand-soft);padding:12px 14px;border-radius:0 6px 6px 0}
table{border-collapse:collapse;width:100%} th,td{text-align:left;padding:6px 9px;border-bottom:1px solid var(--line);font-size:13px;vertical-align:top}
th{background:var(--panel);font-weight:600}
.links{display:flex;flex-wrap:wrap;gap:8px}.links a{border:1px solid var(--line);border-radius:6px;padding:6px 10px;text-decoration:none;color:var(--brand)}
footer{color:var(--muted);font-size:12px;padding:24px 0}
"""

_TEMPLATE = (
    """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{TITLE}}</title>
<style>"""
    + _BASE_CSS
    + """
.masthead{display:flex;gap:16px;align-items:center}.logo{width:64px;height:64px;object-fit:contain}
.meta{color:var(--muted);font-size:12px;word-break:break-all}
nav.tabs{position:sticky;top:0;z-index:5;background:var(--bg);border-bottom:1px solid var(--line)}
nav.tabs .wrap{display:flex;gap:4px;overflow-x:auto}
nav.tabs a{padding:10px 12px;color:var(--muted);text-decoration:none;white-space:nowrap;border-bottom:2px solid transparent}
nav.tabs a:hover{color:var(--ink);border-color:var(--brand)}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0}
.card{border:1px solid var(--line);border-radius:8px;padding:12px;background:var(--panel)}
.card span{display:block;color:var(--muted);font-size:12px}.card strong{display:block;font-size:21px;margin-top:4px}
.finding{border:1px solid var(--line);border-left-width:4px;border-radius:6px;padding:10px 12px;margin:8px 0;background:var(--panel)}
.finding.error{border-left-color:var(--err)}.finding.warning{border-left-color:var(--warn)}.finding.info{border-left-color:var(--info)}.finding.ok{border-left-color:var(--ok)}
.finding .sev{font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.05em}
.finding.error .sev{color:var(--err)}.finding.warning .sev{color:var(--warn)}.finding.info .sev{color:var(--info)}.finding.ok .sev{color:var(--ok)}
.finding details{margin-top:4px;color:var(--muted);font-size:12px}
.controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:8px 0}
.controls label{font-size:12px;color:var(--muted);display:flex;gap:6px;align-items:center}
input,select,button{font:inherit;font-size:13px;color:var(--ink);background:var(--bg);border:1px solid var(--line);border-radius:6px;padding:5px 8px}
button{cursor:pointer}button:hover{border-color:var(--brand)}
#vwrap{position:relative;overflow:auto;height:460px;border:1px solid var(--line);border-radius:6px;background:var(--bg)}
#vcanvas{position:absolute;left:0;top:0}
.tip{position:fixed;pointer-events:none;background:var(--ink);color:var(--bg);font-size:12px;padding:6px 8px;border-radius:5px;display:none;z-index:20;max-width:320px}
.chart{width:100%;height:220px;border:1px solid var(--line);border-radius:6px;background:var(--bg)}
.tablebox{overflow-x:auto;border:1px solid var(--line);border-radius:6px}
.tablebox th{cursor:pointer;user-select:none;white-space:nowrap;position:sticky;top:0}
.tablebox td{white-space:nowrap;font-variant-numeric:tabular-nums}
.pager{display:flex;gap:8px;align-items:center;font-size:12px;color:var(--muted);margin-top:6px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:16px}
figure{margin:0;border:1px solid var(--line);border-radius:8px;padding:10px;background:#fff}
figure img{max-width:100%;height:auto;display:block}figcaption{color:#5b6773;font-size:12px;margin-top:6px}
.muted{color:var(--muted)} code{font-size:12px}
canvas.heat{border:1px solid var(--line);border-radius:6px;max-width:100%}
@media (max-width:640px){h1{font-size:21px}#vwrap{height:360px}}
</style>
</head>
<body>
<header class="top"><div class="wrap masthead">{{LOGO}}
<div><h1>{{TITLE}}</h1><div class="meta">Input: {{INPUT}} &middot; Generated {{GENERATED}} &middot; MSATK {{VERSION}}</div></div>
<div style="margin-left:auto"><button id="themeBtn" title="Toggle light/dark">◐</button></div></div></header>
<nav class="tabs"><div class="wrap">
<a href="#overview">Overview</a><a href="#qc">QC</a><a href="#viewer">Alignment</a><a href="#sites">Sites</a>
<a href="#sequences">Sequences</a><a href="#pairwise">Pairwise</a><a href="#figures">Figures</a><a href="#repro">Reproducibility</a>
</div></nav>
<main class="wrap">
<section id="overview"><h2>Overview</h2>
<div class="cards">{{CARDS}}</div>
<p class="callout">{{INTERPRETATION}}</p>
<h3>Recommended next steps</h3><ul>{{STEPS}}</ul>
<h3>All summary statistics</h3><div id="summaryTable"></div>
</section>
<section id="qc"><h2>Quality assessment</h2><p class="muted">Findings are ranked by severity. Expand a finding to see affected sequences or sites.</p><div id="findings"></div></section>
<section id="viewer"><h2>Alignment viewer</h2>
<div id="viewerNote" class="muted"></div>
<div class="controls">
<label>Search IDs <input id="vsearch" placeholder="substring" size="14"></label>
<label>Go to site <input id="vgoto" type="number" min="1" style="width:90px"></label>
<label>Zoom <input id="vzoom" type="range" min="3" max="22" value="14"></label>
<label>Colour <select id="vmode"><option value="residue">Residue</option><option value="diff">Differences from consensus</option><option value="gaps">Gaps / missing only</option></select></label>
<label><input id="vletters" type="checkbox" checked> Letters</label>
</div>
<div id="vwrap"><div id="vspacer"></div><canvas id="vcanvas"></canvas></div>
</section>
<section id="sites"><h2>Site statistics</h2>
<div class="controls"><label>Track <select id="trackSel"></select></label><span class="muted">Click the track to jump the viewer to that site.</span></div>
<svg id="track" class="chart" preserveAspectRatio="none"></svg>
<h3>Per-site table</h3><div id="siteTable"></div></section>
<section id="sequences"><h2>Sequence statistics</h2><div id="seqTable"></div></section>
<section id="pairwise"><h2>Pairwise identity</h2><div id="pairNote" class="muted"></div><div class="grid2"><div><canvas id="heat" class="heat"></canvas></div><div><svg id="embed" class="chart" style="height:320px;display:none"></svg></div></div></section>
<section id="figures"><h2>Figures</h2><div id="figs" class="grid2"></div></section>
<section id="repro"><h2>Reproducibility</h2>
<div id="provTable"></div>
{{LINKS}}
<h3>Metric dictionary</h3><p class="muted">Definition, formula, missing-data policy and applicability of every reported statistic.</p><div id="metricTable"></div>
</section>
<footer>Generated by MSATK {{VERSION}}. This report is fully offline: it makes no network requests.</footer>
</main>
<div class="tip" id="tip"></div>
<script id="msatk-data" type="application/json">{{DATA}}</script>
<script>
(function(){
"use strict";
const D=JSON.parse(document.getElementById('msatk-data').textContent);
const $=id=>document.getElementById(id);
const tip=$('tip');
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
function fmt(v,k){if(v===null||v===undefined)return'NA';if(typeof v==='boolean')return v?'yes':'no';
 if(typeof v==='number'){if(Number.isInteger(v))return v.toLocaleString();
  if(k&&/fraction|identity|occupancy|gc_content|support/.test(k)&&!/skew|_z$/.test(k)&&v>=0&&v<=1)return (v*100).toFixed(1)+'%';
  return Math.abs(v)>=1e4||Math.abs(v)<1e-3&&v!==0?v.toExponential(3):(+v.toPrecision(4)).toString();}
 if(Array.isArray(v))return v.join('; ');return String(v);}
function esc(s){return String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
function showTip(e,html){tip.innerHTML=html;tip.style.display='block';const x=Math.min(e.clientX+14,innerWidth-330);tip.style.left=x+'px';tip.style.top=(e.clientY+14)+'px';}
function hideTip(){tip.style.display='none';}
// theme toggle
$('themeBtn').onclick=()=>{const r=document.documentElement;const dark=r.dataset.theme?r.dataset.theme==='dark':matchMedia('(prefers-color-scheme: dark)').matches;r.dataset.theme=dark?'light':'dark';try{localStorage.setItem('msatk-theme',r.dataset.theme)}catch(e){}drawAll();};
try{const t=localStorage.getItem('msatk-theme');if(t)document.documentElement.dataset.theme=t;}catch(e){}

// ---------- generic table
function table(el,cols,rows,opts){opts=opts||{};let sortCol=-1,asc=true,page=0,q='';const per=opts.per||50;
 el.innerHTML='<div class="controls"><input placeholder="Filter rows" size="20"><button>Download CSV</button><span class="muted"></span></div><div class="tablebox"><table><thead></thead><tbody></tbody></table></div><div class="pager"><button>‹ Prev</button><span></span><button>Next ›</button></div>';
 const [inp,btn,count]=el.querySelectorAll('.controls > *');const thead=el.querySelector('thead'),tbody=el.querySelector('tbody');const [prev,info,next]=el.querySelector('.pager').children;
 thead.innerHTML='<tr>'+cols.map((c,i)=>'<th data-i="'+i+'">'+esc(opts.labels?opts.labels[i]:c)+'</th>').join('')+'</tr>';
 function filtered(){let r=rows;if(q){const ql=q.toLowerCase();r=r.filter(x=>x.some(v=>v!==null&&String(v).toLowerCase().includes(ql)));}
  if(sortCol>=0){r=r.slice().sort((a,b)=>{const x=a[sortCol],y=b[sortCol];if(x===y)return 0;if(x===null)return 1;if(y===null)return -1;return (x<y?-1:1)*(asc?1:-1);});}return r;}
 function render(){const r=filtered();const pages=Math.max(1,Math.ceil(r.length/per));page=Math.min(page,pages-1);
  tbody.innerHTML=r.slice(page*per,page*per+per).map(x=>'<tr>'+x.map((v,i)=>'<td>'+esc(opts.raw?(v===null?'':v):fmt(v,cols[i]))+'</td>').join('')+'</tr>').join('');
  info.textContent='Page '+(page+1)+' of '+pages;count.textContent=r.length.toLocaleString()+' rows';prev.disabled=page===0;next.disabled=page>=pages-1;}
 thead.onclick=e=>{const i=+e.target.dataset.i;if(isNaN(i))return;asc=sortCol===i?!asc:true;sortCol=i;render();};
 inp.oninput=()=>{q=inp.value;page=0;render();};prev.onclick=()=>{page--;render();};next.onclick=()=>{page++;render();};
 btn.onclick=()=>{const r=filtered();const lines=[cols.join(',')].concat(r.map(x=>x.map(v=>v===null?'':/[",\\n]/.test(String(v))?'"'+String(v).replace(/"/g,'""')+'"':v).join(',')));
  const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([lines.join('\\n')],{type:'text/csv'}));a.download=(opts.name||'table')+'.csv';a.click();};
 render();}

// ---------- overview / provenance
const prov=['msatk_version','output_schema_version','input','input_sha256','generated_at','python_version','platform','coordinate_convention','pairwise_method','pairwise_sample_seed','detected_format','molecule_type'];
const S=D.summary;
table($('summaryTable'),['statistic','value'],Object.keys(S).filter(k=>!prov.includes(k)&&k!=='interpretation'&&k!=='recommended_next_steps').map(k=>[k,fmt(S[k],k)]),{raw:true,name:'summary',per:25});
table($('provTable'),['field','value'],prov.filter(k=>k in S).map(k=>[k,S[k]===null?'NA':String(S[k])]),{raw:true,name:'provenance',per:20});
table($('metricTable'),['name','level','description','formula','value_range','undefined','missing_policy'],D.metrics.map(m=>[m.name,m.level,m.description,m.formula,m.value_range,m.undefined,m.missing_policy]),{raw:true,name:'metric_dictionary',per:15});

// ---------- findings
$('findings').innerHTML=D.findings.map(f=>'<div class="finding '+f.severity+'"><div class="sev">'+esc(f.severity)+' · '+esc(f.code)+'</div><div>'+esc(f.message)+'</div>'+(f.suggestion?'<div class="muted">'+esc(f.suggestion)+'</div>':'')+
 (f.items&&f.items.length?'<details><summary>'+f.items.length+' affected '+esc(f.scope)+'(s)</summary>'+esc(f.items.slice(0,300).join(', '))+(f.items.length>300?' …':'')+'</details>':'')+'</div>').join('');

// ---------- site tables & track
const SC=D.sites.columns,SR=D.sites.rows;
table($('siteTable'),SC,SR,{name:'per_site_stats'});
table($('seqTable'),D.sequences.columns,D.sequences.rows,{name:'per_sequence_stats'});
const trackKeys=['entropy','normalized_entropy','conservation_score','consensus_support','gap_fraction','occupancy','gc_content'].filter(k=>SC.includes(k));
$('trackSel').innerHTML=trackKeys.map(k=>'<option>'+k+'</option>').join('');
function drawTrack(){const svg=$('track');const k=$('trackSel').value;const ci=SC.indexOf(k);const W=svg.clientWidth||900,H=svg.clientHeight||220,pad=32;
 const vals=SR.map(r=>r[ci]);const fin=vals.filter(v=>v!==null);const max=Math.max(1e-9,...fin,/fraction|occupancy|support|score|normalized|gc/.test(k)?1:0);
 svg.setAttribute('viewBox','0 0 '+W+' '+H);const n=vals.length;const x=i=>pad+(W-pad-8)*(n<=1?0.5:i/(n-1));const y=v=>H-20-(H-34)*(v/max);
 let d='',pen=false;vals.forEach((v,i)=>{if(v===null){pen=false;return;}d+=(pen?'L':'M')+x(i).toFixed(1)+','+y(v).toFixed(1);pen=true;});
 const ink=css('--muted'),brand=k==='gap_fraction'?css('--warn'):css('--brand');
 svg.innerHTML='<line x1="'+pad+'" y1="'+(H-20)+'" x2="'+(W-8)+'" y2="'+(H-20)+'" stroke="'+ink+'" stroke-width="0.6"/>'+
 '<text x="4" y="14" font-size="10" fill="'+ink+'">'+fmt(max,k)+'</text><text x="4" y="'+(H-22)+'" font-size="10" fill="'+ink+'">0</text>'+
 '<text x="'+pad+'" y="'+(H-5)+'" font-size="10" fill="'+ink+'">1</text><text x="'+(W-8)+'" y="'+(H-5)+'" font-size="10" text-anchor="end" fill="'+ink+'">'+n+'</text>'+
 '<path d="'+d+'" fill="none" stroke="'+brand+'" stroke-width="'+(n>800?0.8:1.5)+'"/><line id="cursor" y1="8" y2="'+(H-20)+'" stroke="'+ink+'" stroke-dasharray="3,3" visibility="hidden"/>';
 const idx=e=>{const r=svg.getBoundingClientRect();const fx=(e.clientX-r.left)/r.width*W;return Math.max(0,Math.min(n-1,Math.round((fx-pad)/((W-pad-8)/Math.max(1,n-1)))));};
 svg.onmousemove=e=>{const i=idx(e);const c=svg.querySelector('#cursor');c.setAttribute('x1',x(i));c.setAttribute('x2',x(i));c.setAttribute('visibility','visible');
  showTip(e,'Site '+SR[i][0]+'<br>'+k+': '+fmt(vals[i],k)+'<br>class: '+SR[i][SC.indexOf('site_class')]);};
 svg.onmouseleave=hideTip;svg.onclick=e=>{gotoSite(SR[idx(e)][0]);location.hash='#viewer';};}
$('trackSel').onchange=drawTrack;

// ---------- alignment viewer
const V=D.viewer;let vrows=[],cw=14,ch=16;const labelW=150,headH=34;
const vwrap=$('vwrap'),vc=$('vcanvas'),vs=$('vspacer');
function hexLum(h){const n=parseInt(h.slice(1),16);return (0.299*(n>>16)+0.587*(n>>8&255)+0.114*(n&255))/255;}
function cellColor(c,j,mode){const T=D.theme;if(c==='-'||c==='.')return css('--panel');if(c==='?')return '#ced4da';
 if(mode==='gaps')return css('--line');if(mode==='diff'&&c===V.consensus[j])return css('--bg');return T[c]||'#868e96';}
function setupViewer(){if(!V){$('viewerNote').textContent='The interactive viewer is disabled for alignments larger than '+D.viewerLimit.toLocaleString()+' cells; use `msatk plot` to render a region.';vwrap.style.display='none';return;}
 const q=$('vsearch').value.toLowerCase();vrows=V.ids.map((_,i)=>i).filter(i=>!q||V.ids[i].toLowerCase().includes(q));
 cw=+$('vzoom').value;ch=Math.max(cw+2,8);vs.style.width=(labelW+V.consensus.length*cw)+'px';vs.style.height=(headH+(vrows.length+1)*ch+ch)+'px';
 $('viewerNote').textContent=vrows.length+' of '+V.ids.length+' sequences shown · '+V.consensus.length.toLocaleString()+' sites · top row is the majority consensus';drawViewer();}
function drawViewer(){if(!V)return;const dpr=window.devicePixelRatio||1;const W=vwrap.clientWidth,H=vwrap.clientHeight;vc.width=W*dpr;vc.height=H*dpr;vc.style.width=W+'px';vc.style.height=H+'px';
 vc.style.left=vwrap.scrollLeft+'px';vc.style.top=vwrap.scrollTop+'px';const g=vc.getContext('2d');g.setTransform(dpr,0,0,dpr,0,0);
 g.fillStyle=css('--bg');g.fillRect(0,0,W,H);const mode=$('vmode').value,letters=$('vletters').checked&&cw>=8;
 const L=V.consensus.length;const c0=Math.floor(vwrap.scrollLeft/cw),c1=Math.min(L,c0+Math.ceil((W-labelW)/cw)+1);
 const r0=Math.floor(vwrap.scrollTop/ch),nVis=Math.ceil((H-headH-ch)/ch)+1;const offX=labelW-(vwrap.scrollLeft%cw),offY=headH+ch-(vwrap.scrollTop%ch);
 g.font=(Math.max(7,cw-3))+'px ui-monospace,Menlo,Consolas,monospace';g.textAlign='center';g.textBaseline='middle';
 const drawRow=(seq,y)=>{for(let j=c0;j<c1;j++){const c=seq[j];const x=offX+(j-c0)*cw;const col=cellColor(c,j,mode);g.fillStyle=col;g.fillRect(x,y,cw-(cw>6?1:0),ch-1);
  if(letters&&c!=='-'&&c!=='.'){g.fillStyle=(col.startsWith('#')&&hexLum(col)<0.45)?'#fff':'#111';if(mode==='diff'&&c===V.consensus[j])g.fillStyle=css('--muted');g.fillText(c,x+cw/2,y+ch/2);}}};
 g.save();g.beginPath();g.rect(labelW,headH+ch,W-labelW,H-headH-ch);g.clip();
 for(let k=0;k<nVis&&r0+k<vrows.length;k++)drawRow(V.rows[vrows[r0+k]],offY+k*ch);g.restore();
 g.save();g.beginPath();g.rect(labelW,0,W-labelW,headH+ch);g.clip();
 const step=[1,5,10,20,50,100,200,500,1000,5000].find(s=>s*cw>=48)||10000;g.fillStyle=css('--muted');g.font='10px sans-serif';
 for(let j=c0;j<c1;j++){if((j+1)%step===0||j===0){const x=offX+(j-c0)*cw+cw/2;g.fillText(String(j+1),x,12);g.fillRect(x,20,1,6);}}
 g.font=(Math.max(7,cw-3))+'px ui-monospace,Menlo,Consolas,monospace';drawRow(V.consensus,headH);g.restore();
 g.fillStyle=css('--panel');g.fillRect(0,0,labelW,H);g.fillStyle=css('--ink');g.textAlign='right';g.font='11px sans-serif';
 g.fillText('consensus',labelW-6,headH+ch/2);g.save();g.beginPath();g.rect(0,headH+ch,labelW,H);g.clip();
 for(let k=0;k<nVis&&r0+k<vrows.length;k++){let id=V.ids[vrows[r0+k]];if(id.length>22)id=id.slice(0,21)+'…';g.fillText(id,labelW-6,offY+k*ch+ch/2);}g.restore();}
function gotoSite(s){if(!V)return;vwrap.scrollLeft=Math.max(0,(s-1)*cw-(vwrap.clientWidth-labelW)/2);drawViewer();}
if(V){vwrap.addEventListener('scroll',()=>requestAnimationFrame(drawViewer));
 vwrap.addEventListener('mousemove',e=>{const r=vwrap.getBoundingClientRect();const x=e.clientX-r.left,y=e.clientY-r.top;if(x<labelW||y<headH){hideTip();return;}
  const j=Math.floor((x-labelW+vwrap.scrollLeft)/cw);if(j>=V.consensus.length){hideTip();return;}
  const site=SR[j];const info='Site '+(j+1)+' · consensus '+V.consensus[j]+'<br>entropy '+fmt(site[SC.indexOf('entropy')])+' · gaps '+fmt(site[SC.indexOf('gap_fraction')],'gap_fraction');
  if(y<headH+ch){showTip(e,info);return;}const ri=Math.floor((y-headH-ch+vwrap.scrollTop)/ch);if(ri>=vrows.length){hideTip();return;}
  const i=vrows[ri];showTip(e,'<b>'+esc(V.ids[i])+'</b><br>residue '+esc(V.rows[i][j])+'<br>'+info);});
 vwrap.addEventListener('mouseleave',hideTip);
 ['vsearch','vzoom','vmode','vletters'].forEach(id=>$(id).addEventListener('input',setupViewer));
 $('vgoto').addEventListener('change',e=>gotoSite(+e.target.value));}

// ---------- pairwise heatmap & embedding
const ramp=['#440154','#3b528b','#21918c','#5ec962','#fde725'];
function rampColor(v){if(v===null)return css('--panel');const t=Math.max(0,Math.min(1,v))*(ramp.length-1);const i=Math.min(ramp.length-2,Math.floor(t)),f=t-i;
 const a=parseInt(ramp[i].slice(1),16),b=parseInt(ramp[i+1].slice(1),16);const m=s=>Math.round(((a>>s)&255)*(1-f)+((b>>s)&255)*f);return 'rgb('+m(16)+','+m(8)+','+m(0)+')';}
function drawHeat(){const H=D.heatmap,P=D.pairSummary,c=$('heat');
 if(!P){$('pairNote').textContent='Pairwise comparison was not computed (pairwise method: none, or fewer than two sequences).';c.style.display='none';return;}
 $('pairNote').innerHTML='Method: <b>'+esc(P.pairwise_method)+'</b> · '+P.pairwise_pairs_evaluated.toLocaleString()+' of '+P.pairwise_pairs_total.toLocaleString()+' pairs · mean identity '+fmt(P.mean_pairwise_identity,'identity')+' · pairs without comparable sites: '+P.pairwise_pairs_without_comparable_sites+(H?'':' · heatmap shown only for exact mode with ≤ 400 sequences.');
 if(!H){c.style.display='none';return;}const n=H.ids.length,px=Math.max(1,Math.floor(Math.min(520,c.parentElement.clientWidth||520)/n));c.width=c.height=n*px;
 const g=c.getContext('2d');for(let i=0;i<n;i++)for(let j=0;j<n;j++){g.fillStyle=rampColor(i===j?1:H.values[i][j]);g.fillRect(j*px,i*px,px,px);}
 c.onmousemove=e=>{const r=c.getBoundingClientRect();const j=Math.floor((e.clientX-r.left)/r.width*n),i=Math.floor((e.clientY-r.top)/r.height*n);
  if(i<0||j<0||i>=n||j>=n)return;showTip(e,esc(H.ids[i])+' × '+esc(H.ids[j])+'<br>identity '+fmt(i===j?1:H.values[i][j],'identity'));};c.onmouseleave=hideTip;}
function drawEmbed(){const E=D.embeddings;const svg=$('embed');if(!E.length)return;svg.style.display='block';const W=svg.clientWidth||400,H=320,p=28;
 const xs=E.map(e=>e.embedding_1),ys=E.map(e=>e.embedding_2);const [x0,x1,y0,y1]=[Math.min(...xs),Math.max(...xs),Math.min(...ys),Math.max(...ys)];
 const X=v=>p+(W-2*p)*(x1>x0?(v-x0)/(x1-x0):0.5),Y=v=>H-p-(H-2*p)*(y1>y0?(v-y0)/(y1-y0):0.5);svg.setAttribute('viewBox','0 0 '+W+' '+H);
 svg.innerHTML='<text x="'+p+'" y="16" font-size="11" fill="'+css('--muted')+'">'+esc(String(E[0].method).toUpperCase())+' embedding</text>'+E.map(e=>'<circle cx="'+X(e.embedding_1)+'" cy="'+Y(e.embedding_2)+'" r="3.5" fill="'+css('--brand')+'" fill-opacity="0.75"><title>'+esc(e.sequence_id)+'</title></circle>').join('');}

// ---------- figures
$('figs').innerHTML=D.figures.length?D.figures.map(f=>'<figure><img loading="lazy" src="'+f.src+'" alt="'+esc(f.name)+'"><figcaption>'+esc(f.name)+'</figcaption></figure>').join(''):'<p class="muted">No static figures in this report (install the plots extra or use --plots png).</p>';

function drawAll(){drawTrack();drawViewer();drawHeat();drawEmbed();}
setupViewer();drawTrack();drawHeat();drawEmbed();
let rt;addEventListener('resize',()=>{clearTimeout(rt);rt=setTimeout(drawAll,120);});
})();
</script>
</body>
</html>
"""
)
