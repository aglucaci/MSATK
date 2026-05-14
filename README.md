# MSATK

![Tests](https://github.com/yourname/msatk/actions/workflows/tests.yml/badge.svg)
![Docs](https://github.com/yourname/msatk/actions/workflows/docs.yml/badge.svg)
![PyPI](https://img.shields.io/pypi/v/msatk)
![Python](https://img.shields.io/pypi/pyversions/msatk)
![License](https://img.shields.io/github/license/yourname/msatk)

**MSATK** is a Python and command-line toolkit for profiling multiple sequence alignments. It converts nucleotide, protein, and codon-aware alignments into quality-control summaries, per-sequence and per-site statistics, publication-ready plots, codon usage metrics, RSCU tables, and standalone HTML reports.

> **MSATK: one command to profile, visualize, and report multiple sequence alignments.**

## Quick Start

```bash
pip install msatk
msatk profile alignment.fasta
```

MSATK automatically detects the input format and molecule type, computes summary statistics, generates plots when plotting dependencies are available, and writes a complete report.

```text
alignment_msatk/
├── report.html
├── summary.json
├── parameters.yaml
├── msatk.log
├── qc_warnings.txt
├── tables/
└── figures/
```

Try the bundled demo:

```bash
msatk demo
```

## Why MSATK?

- One-command alignment profiling
- DNA, RNA, protein, and codon-aware modes
- Publication-ready PNG plots with optional plotting dependencies
- CSV, JSON, HTML, and Markdown outputs
- Codon usage, RSCU, GC1/GC2/GC3, and stop codon summaries
- Protein amino-acid, residue-class, and hydrophobicity summaries
- QC warnings with plain-language interpretation
- Python API for notebooks and workflows
- Stable output files for Snakemake, Nextflow, and HPC pipelines

## Install Extras

```bash
pip install "msatk[all]"
```

Useful extras:

```bash
pip install "msatk[dataframes]"
pip install "msatk[plots]"
pip install "msatk[embed]"
```

Developer install:

```bash
pip install -e ".[dev,all]"
```

## CLI

```bash
msatk profile alignment.fasta
msatk qc alignment.fasta
msatk codon cds_alignment.fasta
msatk protein protein_alignment.faa
msatk embed alignment.fasta --method pca
msatk report alignment_msatk/
msatk demo
```

Example terminal output:

```text
MSATK alignment profile complete

Input: alignment.fasta
Detected type: codon
Sequences: 248
Alignment length: 3,642
Gap fraction: 4.8%
Variable sites: 712
Mean pairwise identity: 91.4%

Results written to:
alignment_msatk/
Open report:
alignment_msatk/report.html
```

## Python API

```python
from msatk import AlignmentProfiler, CodonProfiler, ProteinProfiler, profile_alignment

profiler = AlignmentProfiler("alignment.fasta")
summary = profiler.summary()
per_sequence = profiler.per_sequence_stats()
per_site = profiler.per_site_stats()

profiler.plot_gap_profile()
profiler.plot_entropy()
profiler.write_report("report.html")

results = profile_alignment("alignment.fasta", outdir="alignment_msatk")

codon = CodonProfiler("cds_alignment.fasta")
rscu = codon.rscu()

protein = ProteinProfiler("protein_alignment.faa")
aa = protein.amino_acid_composition()
```

When pandas is installed, table-like API methods return pandas DataFrames. In minimal environments, MSATK falls back to lists of dictionaries.

## Roadmap

- `0.1`: one-command profiling, summary statistics, QC, tables, report
- `0.2`: hardened format support, better terminal output, example datasets
- `0.3`: codon-aware release with RSCU, GC1/GC2/GC3, stop codon detection
- `0.4`: visualization themes, SVG output, richer report customization
- `0.5`: embeddings, clustering, outlier detection
- `1.0`: stable API/output schema, PyPI/Bioconda, containers, workflow examples
