---
title: "MSATK: a Python toolkit for statistical, visual, and codon-aware profiling of multiple sequence alignments"
tags:
  - Python
  - bioinformatics
  - multiple sequence alignment
  - phylogenetics
  - codon usage
  - quality control
authors:
  - name: Alexander Lucaci
    affiliation: 1
affiliations:
  - name: MSATK Project
    index: 1
date: 14 May 2026
bibliography: references.bib
---

# Abstract

Multiple sequence alignments are central to molecular evolution, comparative genomics, protein analysis, phylogenetics, and genomic surveillance. However, researchers often need to combine ad hoc scripts and several specialized tools to summarize alignment quality, composition, conservation, codon usage, and sequence-level outliers. We present MSATK, a Python and command-line toolkit for profiling nucleotide, protein, codon-aware, and alignment-derived sequencing inputs. MSATK automatically detects common alignment formats, infers molecule type, computes alignment-wide, per-sequence, and per-site statistics, generates quality-control warnings, and writes machine-readable tables together with standalone HTML and Markdown reports. For coding alignments, MSATK reports codon usage, relative synonymous codon usage, codon-position summaries, stop codon warnings, frame-related warnings, and translated amino-acid composition. Protein alignments are summarized using amino-acid composition, residue-class summaries, hydrophobicity summaries, conservation, and pairwise identity metrics. MSATK is designed for interactive command-line use, Python notebooks, and automated workflows requiring stable output files and reproducible metadata. The package provides a lightweight and extensible framework for alignment quality control and exploratory analysis prior to phylogenetic inference, selection analysis, protein-family analysis, or downstream comparative genomics.

# Introduction

Multiple sequence alignments (MSAs) are a core data structure in modern bioinformatics. They support phylogenetic inference, molecular evolution, comparative genomics, protein family analysis, metagenomics, and viral genomic surveillance. Before downstream analyses, researchers typically need to evaluate alignment completeness, gap structure, sequence quality, site-level conservation, compositional bias, pairwise similarity, and, for coding sequences, codon-level properties. These summaries are often generated through ad hoc scripts or by combining multiple specialized tools, making results difficult to reproduce and compare across projects.

Existing software provides excellent support for alignment construction, visualization, trimming, phylogenetic inference, and molecular evolution. Tools such as MAFFT [@katoh2013mafft], MUSCLE [@edgar2004muscle], and Clustal Omega [@sievers2011clustal] construct alignments; Jalview [@waterhouse2009jalview] and AliView [@larsson2014aliview] support interactive inspection; trimAl [@capella2009trimal] assists with alignment trimming; IQ-TREE [@minh2020iqtree], RAxML [@stamatakis2014raxml], HyPhy [@pond2020hyphy], and PAML [@yang2007paml] support downstream inference. These tools are not intended to provide a single standardized, report-oriented summary of alignment quality, composition, codon usage, and reproducible output metadata. As a result, users frequently bridge the gap between alignment generation and downstream inference using project-specific scripts.

MSATK addresses this gap by providing a "FastQC-like" profiling layer for multiple sequence alignments. Given an input alignment, MSATK produces a structured output directory containing summary metadata, quality-control warnings, tabular statistics, optional plots, and shareable reports. The goal is not to replace alignment algorithms or downstream inference software, but to help users rapidly understand whether an alignment is complete, coherent, codon-aware when expected, and suitable for subsequent analyses.

# Software Overview

MSATK is both a Python package and a command-line application. The primary command is:

```bash
msatk profile alignment.fasta --out msatk_results/
```

The command automatically detects the input format and molecule type, validates alignment structure, computes descriptive statistics, generates warnings, and writes standardized outputs. A bundled demonstration command is also provided:

```bash
msatk demo
```

MSATK can also be used programmatically:

```python
from msatk import MSATK

profiler = MSATK("alignment.fasta")
summary = profiler.summary()
per_sequence = profiler.per_sequence_stats()
per_site = profiler.per_site_stats()
profiler.write_outputs("msatk_results", force=True)
```

For compatibility with earlier design sketches, `AlignmentProfiler` remains available as an alias for `MSATK`.

## Supported Inputs

MSATK supports common multiple-alignment and alignment-derived formats:

- FASTA and aligned FASTA
- PHYLIP and relaxed PHYLIP
- NEXUS matrix blocks
- CLUSTAL
- Stockholm
- A3M
- MAF
- SAM
- BAM and CRAM-derived alignment summaries through optional `pysam`

Supported sequence and analysis modes include nucleotide alignments, amino-acid alignments, codon alignments, translated CDS-like inputs, and mixed or auto-detected mode. MSATK automatically reports file format, molecule type, alignment length, whether sequence lengths are consistent, whether the length is divisible by three, whether stop codons occur, whether frame-related warnings are present, and whether the alignment appears codon-aware.

## Outputs

Each profiling run writes a reproducible output directory. The default layout includes:

```text
msatk_results/
|-- report.html
|-- report.md
|-- summary.json
|-- parameters.yaml
|-- msatk.log
|-- qc_warnings.txt
|-- assets/
|-- tables/
`-- figures/
```

Core tables include `alignment_summary.csv`, `per_sequence_stats.csv`, `per_site_stats.csv`, `pairwise_identity_matrix.csv`, `distance_matrix.csv`, `composition_summary.csv`, `entropy_by_site.csv`, `sequence_embeddings.csv`, and `qc_flags.csv`. Codon-aware runs additionally write codon usage, RSCU, codon-position GC, codon-position entropy, stop codon, and translated amino-acid summaries. Protein runs write amino-acid composition, residue-class, and hydrophobicity summaries.

# Methods and Implementation

MSATK is implemented in Python and organized as a modular package with independent layers for input parsing, validation, statistics, codon analysis, protein analysis, plotting, embeddings, report generation, and workflow integration. The core package uses lightweight internal parsers for common MSA formats, allowing basic profiling without requiring large dependencies. Optional dependencies extend the user experience: pandas [@mckinney2010pandas] provides notebook-friendly DataFrame outputs, Matplotlib [@hunter2007matplotlib] enables PNG figures, scikit-learn [@pedregosa2011scikit] enables dimensionality-reduction embeddings, and pysam [@pysam] enables BAM and CRAM parsing through HTSlib-compatible interfaces.

The internal workflow is:

```text
Input reader -> validation -> detection -> statistics -> plots -> reports -> output writer
```

The validation layer supports permissive and strict modes. In permissive mode, MSATK warns and continues where scientifically reasonable, for example by padding unequal sequence lengths for site-level summaries. In strict mode, duplicate identifiers, invalid characters, and unequal sequence lengths raise explicit validation errors. The command-line interface exposes these modes through `--permissive` and `--strict`.

## Statistics Engine

MSATK computes alignment-wide, per-sequence, and per-site summaries. Alignment-level metrics include sequence count, alignment length, gap fraction, missingness, conserved sites, variable sites, parsimony-informative sites, singleton sites, mean pairwise identity, mean pairwise distance, entropy summaries, and GC content for nucleotide-like alignments. Per-sequence outputs include ungapped length, gap fraction, ambiguous-character fraction, GC content, unique residue count, mean identity to other sequences, and a simple outlier score. Per-site outputs include dominant character, dominant-character frequency, gap fraction, missing fraction, Shannon entropy [@shannon1948], conservation score, observed-state count, variable-site status, parsimony-informative status, singleton status, and codon position.

## Codon and Protein Summaries

For codon-aware alignments, MSATK reports codon usage, relative synonymous codon usage (RSCU), codon-position GC content, codon-position entropy, stop codons, and translated amino-acid summaries. RSCU summarizes observed codon use relative to the expected use within synonymous codon families [@sharp1986rscu]. MSATK also reports frame-related warnings, including alignment lengths that are not divisible by three, sequence ungapped lengths that are not divisible by three, and partial gaps inside codon triplets.

For protein alignments, MSATK reports amino-acid composition, residue-class composition, hydrophobicity summaries using the Kyte-Doolittle scale [@kyte1982hydropathy], site-level entropy, conservation, and pairwise identity metrics.

## Reporting and Reproducibility

MSATK reports are generated from structured outputs and include a branded HTML report, Markdown report, CSV tables, JSON metadata, and QC warning files. Each run records the MSATK version, Python version, platform, input path, detected format, detected molecule type, command-line parameters, timestamp, validation mode, and output schema version. This makes the output directory suitable for interactive review and automated workflow consumption.

# Example Application

As a demonstration, consider a small coding-sequence alignment:

```bash
msatk profile tests/data/codon/valid_codon_alignment.fasta --type codon --out codon_results --force
```

MSATK detects the input as a codon-aware nucleotide alignment, reports three sequences and twelve alignment sites, computes low missingness, summarizes variable sites and mean pairwise identity, and writes codon-specific tables including `codon_usage.csv`, `rscu.csv`, `codon_position_gc.csv`, `codon_position_entropy.csv`, `stop_codon_report.csv`, and `translated_amino_acid_summary.csv`. A related fixture containing an internal stop codon triggers a plain-language warning identifying the affected sequence. These outputs provide a compact pre-analysis check before running selection-analysis tools such as HyPhy or PAML.

The same workflow can be applied to protein alignments:

```bash
msatk protein tests/data/protein/protein_with_gaps.faa --out protein_results --force
```

In this mode, MSATK writes amino-acid composition, residue class, hydrophobicity, pairwise identity, entropy, and conservation summaries. For workflow developers, all outputs are deterministic filenames under a single output directory, allowing direct integration with Snakemake, Nextflow, or HPC batch jobs.

# Availability

MSATK is freely available under the MIT License at `https://github.com/yourname/msatk`. The package supports Python 3.10 and later. Source installation and editable development installation are supported with:

```bash
pip install -e ".[dev,all,docs]"
```

MSATK includes a local Conda recipe and a Bioconda submission template. The intended bioinformatics-native installation path is:

```bash
mamba install -c bioconda -c conda-forge msatk
```

Until a Bioconda package is published, users can install Conda-managed dependencies using the provided `environment.yml` and install MSATK from source or PyPI. Documentation, synthetic test datasets, GitHub Actions workflows, and release checklists are included with the repository.

# Documentation and Reproducibility

The repository includes unit tests, integration tests, regression tests, synthetic alignment fixtures, CI workflows, documentation for CLI and Python API usage, Conda packaging files, and example datasets. Test data cover FASTA, aligned FASTA, A3M, PHYLIP, relaxed PHYLIP, CLUSTAL, Stockholm, NEXUS, MAF, SAM, and BAM fixtures, with optional CRAM fixture generation where HTSlib-compatible tooling is available. Reports and tables include versioned metadata and a stable output schema to support reproducible research and workflow automation.

# Discussion

MSATK contributes a standardized profiling and reporting layer for multiple sequence alignments. It lowers the barrier to quality control by turning an alignment into a coherent set of statistics, warnings, tables, figures, and reports with one command. This is useful for teaching, exploratory data analysis, manuscript supplementary materials, and production bioinformatics workflows.

MSATK is not intended to replace alignment construction, alignment visualization, phylogenetic inference, or selection-analysis software. Instead, it sits between alignment generation and downstream inference. By combining gap/missingness summaries, entropy, conservation, composition, codon-aware metrics, protein summaries, pairwise identity, warnings, and reproducible reporting, MSATK helps users identify problematic alignments and generate standardized summaries before committing to computationally expensive or biologically sensitive analyses.

Future development will focus on expanding validated input coverage, improving large-alignment performance, adding richer SVG and interactive plotting options, strengthening embedding and outlier-detection workflows, and providing workflow wrappers for Snakemake, Nextflow, Galaxy, and containerized deployments.

# Conclusion

MSATK provides a lightweight, reproducible, and user-friendly framework for profiling multiple sequence alignments. By combining alignment statistics, quality-control warnings, codon-aware summaries, protein summaries, visualizations, and standalone reports, MSATK helps researchers rapidly assess alignment quality and prepare data for evolutionary, comparative, and functional analyses.

# Acknowledgements

MSATK builds on the broader open-source scientific Python and bioinformatics ecosystem. The author thanks contributors to Python, pandas, Matplotlib, scikit-learn, pysam, and the many alignment and evolutionary-analysis tools that motivate reproducible alignment quality-control workflows.

# Funding

No external funding was received for this work unless otherwise specified by the author.

# Conflict of Interest

The author declares no competing interests.

# References
