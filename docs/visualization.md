# Visualization

## Residue-level alignment figures

```python
from msatk import MSATK

msa = MSATK("proteins.faa")
msa.plot_alignment(
    "figure.svg",
    start=1, end=120,                 # 1-based inclusive region
    order="consensus_identity",       # or "input", "gaps"
    theme="protein",                  # or nucleotide, nucleotide_cvd, grayscale, dict, JSON/YAML file
    tracks=("conservation", "entropy", "gaps"),
    highlight_differences=True,       # fade residues identical to the consensus
    show_consensus=True,
    max_rows=100,
    dpi=300,
)
msa.plot_logo("logo.pdf", start=40, end=80)
```

The output format follows the file suffix (PNG, SVG, PDF, TIFF, JPEG). For large
alignments render a region: whole-genome matrices compressed into one image are
not readable. Letters are drawn automatically when the region is small enough.

Tracks always describe the **whole alignment** at the displayed columns, so they do not
change when rows are subset or truncated. The `conservation` track ("consensus support")
is the fraction of *all* sequences that carry the consensus residue, so a column filled by
only a few sequences cannot look perfectly conserved. `gaps`, `entropy` and `occupancy`
are the per-site gap fraction, normalised entropy and residue occupancy.

Custom colour themes are JSON or YAML mappings of symbol → hex colour:

```yaml
A: "#3fa34d"
C: "#2f6fdf"
G: "#f2a900"
T: "#d7263d"
"-": "#ffffff"
```

Sequence logos show information content (bits) with the Schneider et al. (1986)
small-sample correction (disable with `small_sample_correction=False`).

## Standard figures

`msatk profile` writes per-site gap, entropy and conservation profiles, per-sequence
missingness, composition, pairwise identity/distance distributions, a pairwise
heatmap (exact mode, ≤ 500 sequences), an optional embedding scatter, an alignment
overview and a sequence logo — in every format passed to `--plots`.

## Interactive report

`report.html` is fully offline (no network requests). It contains:

* summary cards, interpretation and severity-ranked QC findings;
* a canvas alignment viewer with scrolling, zoom, ID search, jump-to-site, colouring
  by residue / differences from consensus / gaps, and hover inspection;
* interactive site tracks (click to jump the viewer);
* sortable, filterable, paginated sequence and site tables with CSV download;
* a pairwise identity heatmap and embedding scatter;
* provenance and the full metric dictionary.

`--self-contained` embeds figures as data URIs so the report is a single file.
The viewer is disabled above 3 million cells; use `msatk plot` for those.
