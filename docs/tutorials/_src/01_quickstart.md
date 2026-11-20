# 1. Profile an alignment in one command

**You will:** profile a real protein-family alignment from the command line, read the
summary, look at the QC findings and open the interactive report.

**Data:** the Ras-family seed alignment from Pfam (PF00071; 60 sequences, 229 columns),
bundled with MSATK.

## Get the data

`msatk datasets` lists the bundled tutorial datasets and copies them to a folder.

```bash
msatk datasets ras_pf00071 --out data
```

Compressed input (`.gz`) is read directly, so there is no need to unpack it.

## Check what the file is

`msatk inspect` reports the format, the data model, what the file contains and which
analyses apply, before any computation.

```bash
msatk inspect data/ras_pf00071.sto.gz
```

## Profile it

```bash
msatk profile data/ras_pf00071.sto.gz --out ras_profile --plots png --self-contained
```

The output directory has a fixed structure:

```python
from pathlib import Path

for path in sorted(Path("ras_profile").rglob("*")):
    if path.is_file() and path.parent.name != "assets":
        print(path)
```

## Read the summary

`summary.json` holds the alignment-level statistics and their provenance. Every field is
defined in `metric_dictionary.json`, which is written alongside it.

```python
import json

summary = json.loads(Path("ras_profile/summary.json").read_text())
for key in ["molecule_type", "number_of_sequences", "alignment_length", "gap_fraction",
            "variable_sites", "parsimony_informative_sites", "mean_pairwise_identity",
            "pairwise_method"]:
    print(f"{key:30} {summary[key]}")
```

Statistics that cannot be computed are reported as `null`, never as 0. For example, a
column that contains only gaps is classed as *unobserved* and has no entropy.

## QC findings

QC findings are ranked by severity. Each one says which sequences or sites are affected
and what to do about it.

```python
findings = json.loads(Path("ras_profile/qc_findings.json").read_text())
for f in findings:
    print(f"[{f['severity']:7}] {f['code']:22} n={f['count']:<3} {f['message']}")
    if f["suggestion"]:
        print(f"{'':12}-> {f['suggestion']}")
```

Pfam seed alignments contain insert-state columns that only a few sequences fill, so
`HIGH_GAP_SITES` is expected here.

## The interactive report

Open `ras_profile/report.html` in a browser. It works offline, and with
`--self-contained` it is a single file you can share. It includes:

- a residue viewer with scrolling, zoom, search and difference highlighting;
- site tracks you can click;
- sortable tables you can export as CSV;
- a pairwise identity heatmap, the QC findings, provenance and the metric dictionary.

Static figures are written alongside the report:

<!-- figure: ras_profile/figures/alignment_overview.png | Alignment overview figure -->

## Next

- [Tutorial 2](02_statistics_python.md) works with the same statistics from Python.
- `msatk profile --help` lists the presets. `basic` skips pairwise comparison; `comprehensive` adds exact pairwise identity, embeddings and ML tensors.
