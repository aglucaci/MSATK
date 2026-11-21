# MSATK: Multiple Sequence Alignment Tool-Kit

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="logo/msatk_logo_dark.png">
    <img src="logo/msatk_logo.png" alt="MSATK logo" width="420">
  </picture>
</p>

![Tests](https://github.com/aglucaci/MSATK/actions/workflows/tests.yml/badge.svg)
![Docs](https://github.com/aglucaci/MSATK/actions/workflows/docs.yml/badge.svg)
![PyPI](https://img.shields.io/pypi/v/msatk)
![Python](https://img.shields.io/pypi/pyversions/msatk)
![License](https://img.shields.io/github/license/aglucaci/MSATK)

**MSATK** is a unified framework for **statistical profiling, visualization and
machine-learning representations of multiple sequence alignments**. One command
takes an existing alignment to explicitly defined statistics at alignment, sequence,
site and pair resolution; severity-ranked QC findings; publication-quality figures;
an interactive offline HTML report; and reproducible tensor/image exports for
downstream models.

> **One command, one complete alignment characterization, multiple interoperable outputs, no hidden scientific assumptions.**

### What makes it different

- **Explicit semantics.** Gaps, missing data (`N`/`X`/`?`) and ambiguity codes are
  handled by one molecule-aware policy shared by every metric. Undefined values
  (all-gap columns, non-overlapping pairs) are reported as `null`, never as 0.
  Every output field is defined in a [metric dictionary](docs/metrics.md).
- **Multi-resolution statistics.** 145 defined fields (alignment, sequence, site and pair level, including provenance): completeness, gap
  structure, conservation, entropy/information content, composition and skews,
  consensus- and composition-based outlier scores, pairwise identity, and
  population-genetic statistics (π, Watterson's θ, Tajima's D, haplotype
  diversity, Ts/Tv) validated against scikit-allel.
- **Scales by design.** Vectorised NumPy core; exact pairwise identity via blocked
  matrix products, with a size-aware policy (exact → seeded sampling → off).
- **Interactive offline report.** Scrollable/zoomable residue viewer, clickable
  site tracks, sortable/searchable tables with CSV export, pairwise heatmap,
  provenance — no server, account or network access; optional single-file mode.
- **Publication figures.** Residue-level alignment figures with consensus and
  annotation tracks, sequence logos, custom colour themes; PNG/SVG/PDF/TIFF.
- **ML-ready representations.** Lossless categorical encodings, one-hot with
  explicit ambiguity policies, physicochemical and QC feature channels,
  native, tiled or padded exports with masks and recorded row order (never
  resampled), RGB/label images — each with
  a `manifest.json`, plus NumPy/PyTorch/TensorFlow dataset adapters.

## Quick Start

```bash
mamba install -c bioconda -c conda-forge msatk
msatk profile alignment.fasta
```

MSATK automatically detects the input format and molecule type, computes summary statistics, generates plots when plotting dependencies are available, and writes a complete report.

```text
alignment_msatk/
|-- report.html            interactive offline dashboard
|-- summary.json           statistics + provenance (input SHA-256, versions, parameters)
|-- metric_dictionary.json definitions of every field
|-- qc_findings.json       severity-ranked QC findings
|-- parameters.yaml / run_manifest.json / msatk.log
|-- tables/                per-site, per-sequence, pairwise, gap-block, composition tables
`-- figures/               statistical plots, alignment overview, sequence logo
```

Try the bundled demo:

```bash
msatk demo
```

## Installation

### With Conda

The recommended bioinformatics-native install path is Bioconda:

```bash
mamba install -c bioconda -c conda-forge msatk
```

or:

```bash
conda install -c bioconda -c conda-forge msatk
```

Before the official Bioconda package is available, create an environment with Conda-managed dependencies and install MSATK from PyPI:

```bash
conda env create -f environment.yml
conda activate msatk
```

After Bioconda packaging is live, use:

```bash
conda env create -f environment-bioconda.yml
conda activate msatk
```

### With pip

```bash
pip install msatk
```

### Development install

```bash
git clone https://github.com/aglucaci/MSATK
cd msatk
pip install -e ".[dev,all,docs]"
```

## How To Run Locally

From the repository root:

```bash
pip install -e ".[dev,all,docs]"
```

Run the bundled demo:

```bash
msatk demo --out msatk_demo --force
```

Profile your own alignment:

```bash
msatk profile path/to/alignment.fasta
```

Write to a specific output directory:

```bash
msatk profile path/to/alignment.fasta --out results --force
```

Run codon-aware or protein-specific modes:

```bash
msatk codon path/to/cds_alignment.fasta --out codon_results --force
msatk protein path/to/protein_alignment.faa --out protein_results --force
```

Run tests and checks:

```bash
pytest
ruff check src/msatk tests
ruff format --check src/msatk tests
mkdocs build --strict
```

If you are running without installing the package, set `PYTHONPATH` first:

```bash
# macOS/Linux
PYTHONPATH=src python -m msatk.cli demo --out msatk_demo --force

# Windows PowerShell
$env:PYTHONPATH="src"; python -m msatk.cli demo --out msatk_demo --force
```

## Why MSATK?

- One-command alignment profiling
- DNA, RNA, protein, and codon-aware modes
- Common alignment inputs: FASTA/aligned FASTA, A3M, PHYLIP/relaxed PHYLIP, CLUSTAL, Stockholm, NEXUS, MAF, SAM, BAM, and CRAM
- Publication-ready PNG/SVG/PDF/TIFF figures, residue-level alignment views and logos
- CSV, JSON, HTML, Markdown, NPZ and image outputs
- Codon usage, RSCU, GC1/GC2/GC3, and stop codon summaries
- Protein amino-acid, residue-class, and hydrophobicity summaries
- QC warnings with plain-language interpretation
- Python API for notebooks and workflows
- Stable output files for Snakemake, Nextflow, and HPC pipelines

## Supported Input Types

MSATK distinguishes the alignment **type** (data model), the **method** and the **file
format**, and dispatches each input to analyses that are valid for it:

| Data model | Formats | Analyses |
|---|---|---|
| Sequence alignment (pairwise or multiple) | FASTA, PHYLIP, CLUSTAL, NEXUS, Stockholm, A2M, A3M, MSF, MEGA, SELEX, MASE, EMBOSS | Full statistics, QC, figures, report, ML encodings; pairwise, BLOSUM62, RNA structure, profile-state and partition analyses when applicable |
| Whole-genome blocks | MAF, XMFA, AXT | Block sizes, per-genome coverage, strand and collinearity breakpoints |
| Mapped reads | SAM, BAM, CRAM (pysam) | Mapping rate, MAPQ, CIGAR, quality, depth and breadth |
| Mappings and hits | PAF, PSL, BLAST tabular, GAF | Identity, length, MAPQ/e-value, query and target coverage |
| Sequence graphs | GFA 1, rGFA | Segments, links, components, paths, core/shell/private segments |
| Profile HMMs | HMMER3, HHM | Model metadata, match-state emission entropy |

Gzip-compressed input is accepted. All NCBI genetic codes are available with
`--genetic-code`. Run `msatk inspect file` to see what a file contains and what can be
computed, and `msatk formats` for the full matrix with support levels
([docs/compatibility.md](docs/compatibility.md), generated from the capability registry).
PDB/mmCIF, GBZ, ODGI, HAL and chain files are recognised but not yet supported.

## Install Extras

```bash
pip install "msatk[all]"
```

Useful extras:

```bash
pip install "msatk[dataframes]"
pip install "msatk[plots]"
pip install "msatk[embed]"      # UMAP / t-SNE (PCA needs only NumPy)
pip install "msatk[ngs]"
```

`SAM` is supported without extra dependencies. `BAM` and `CRAM` require `pysam`; install with `pip install "msatk[ngs]"` or `mamba install -c bioconda pysam`.

Developer install:

```bash
pip install -e ".[dev,all]"
```

## Conda And Bioconda Packaging

MSATK is intended to be distributed through PyPI first, then Bioconda:

```text
PyPI first -> local Conda recipe -> Bioconda PR -> CI-tested Conda builds
```

The local Conda recipe is in:

```text
conda-recipe/meta.yaml
```

Build and test it locally:

```bash
mamba create -n conda-build-env -c conda-forge conda-build boa anaconda-client
conda activate conda-build-env
conda build conda-recipe

mamba create -n test-msatk --use-local msatk
conda activate test-msatk
msatk --help
msatk profile --help
```

For Bioconda submission after a PyPI release, use:

```text
conda-recipe/meta.bioconda.yaml
```

Update the PyPI source SHA256, copy it to `recipes/msatk/meta.yaml` in a fork of `bioconda/bioconda-recipes`, and open a pull request.

## CLI

```bash
msatk profile alignment.fasta --out results/                 # full profile + report
msatk profile alignment.fasta --preset comprehensive --plots png,svg,pdf --self-contained
msatk stats alignment.fasta --level site > sites.csv          # summary|site|sequence|pairwise|composition|gaps
msatk qc alignment.fasta
msatk plot alignment.fasta --start 1 --end 300 --out region.svg
msatk logo alignment.fasta --start 1 --end 60 --out logo.pdf
msatk encode alignment.fasta --strategy tiles --height 128 --width 256 --out ml/
msatk image alignment.fasta --kind rgb --scale 4 --out alignment.png
msatk embed alignment.fasta --method pca
msatk codon cds_alignment.fasta
msatk protein protein_alignment.faa
msatk report alignment_msatk/
msatk metrics
msatk demo
```

Expensive steps are opt-in or size-aware: `--pairwise auto|exact|sampled|none`
(default `auto`: exact up to 1,000 sequences, deterministic sampling above) and
`--embeddings none|pca|umap|tsne` (default off). See [docs/CLI.md](docs/CLI.md).

Example terminal output:

```text
MSATK alignment profile complete

Input: alignment.fasta
Detected type: codon
Sequences: 248
Alignment length: 3,642
Gap fraction: 4.8%
Variable sites: 712
Mean pairwise identity: 91.4% (exact)

Results written to:
alignment_msatk/
Open report:
alignment_msatk/report.html
```

## Python API

```python
from msatk import AlignmentProfiler

msa = AlignmentProfiler("alignment.fasta")

summary = msa.summary()          # alignment-level statistics
sites = msa.sites()              # one record per column
sequences = msa.sequences()      # one record per sequence
pairs = msa.pairwise()           # one record per pair
qc = msa.qc()                    # severity-ranked findings

msa.plot_alignment("alignment.svg", start=1, end=300, theme="nucleotide")
msa.plot_logo("logo.pdf", start=1, end=60)

tensor = msa.to_tensor("onehot")             # .array is a NumPy (N, L, C) array
msa.to_image("alignment.png", kind="label")  # lossless label image + manifest
msa.export_encoding("ml/", encoding="onehot", strategy="tiles", height=128, width=256)

msa.report("report.html")                    # interactive single-file report
```

Codon- and protein-specific profilers (`CodonProfiler`, `ProteinProfiler`) add
codon usage, RSCU, GC1/GC2/GC3, stop-codon reports, amino-acid composition,
residue classes and hydropathy. When pandas is installed, tables are returned as
DataFrames; otherwise as lists of dictionaries. See [docs/python_api.md](docs/python_api.md).

## Tutorials

Seven hands-on tutorials in [`docs/tutorials/`](docs/tutorials/index.md) run on data
bundled with the package (`msatk datasets`). The outputs shown were produced by running
the code:

1. [Profile an alignment in one command](docs/tutorials/01_quickstart.md)
2. [Alignment statistics in Python](docs/tutorials/02_statistics_python.md)
3. [QC of a viral coding alignment](docs/tutorials/03_viral_codon_qc.md): 400 influenza HA sequences
4. [Protein families and publication figures](docs/tutorials/04_protein_family_figures.md)
5. [Machine-learning representations](docs/tutorials/05_ml_encodings.md)
6. [RNA structure, profiles, partitions and pairwise alignments](docs/tutorials/06_annotated_alignments.md)
7. [Beyond MSAs: genome blocks, reads, hits, graphs, profile HMMs](docs/tutorials/07_beyond_msas.md)

## Documentation

- [Statistical conventions](docs/statistics.md) and [metric dictionary](docs/metrics.md)
- [Visualization](docs/visualization.md) and [ML encodings](docs/ml_encodings.md)
- [QC findings](docs/qc_warnings.md), [outputs](docs/outputs.md), [benchmarks](docs/benchmarks.md)

## Roadmap

See [docs/ROADMAP.md](docs/ROADMAP.md). Version 0.2 implements the publication-critical
core (statistical semantics, optional/size-aware pairwise, offline interactive
report, alignment figures, ML encodings, metric dictionary, provenance); 1.0 targets
a stable API/schema, external benchmark comparisons and an Application Note.
