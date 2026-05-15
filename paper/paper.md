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

Multiple sequence alignments are central to phylogenetics, molecular evolution, comparative genomics, protein-family analysis, and genomic surveillance. However, researchers often lack a standardized, reproducible way to summarize alignment quality, gap structure, conservation, codon usage, and sequence-level outliers before downstream analysis. We present MSATK, a Python and command-line toolkit for statistical, visual, and codon-aware profiling of multiple sequence alignments. MSATK supports common alignment formats, automatically infers molecule type, computes alignment-wide, per-sequence, and per-site metrics, generates quality-control warnings, and writes machine-readable tables with standalone HTML and Markdown reports. For coding alignments, MSATK reports codon usage, relative synonymous codon usage, codon-position summaries, stop-codon warnings, frame-related warnings, and translated amino-acid composition. Protein alignments are summarized using amino-acid composition, residue-class summaries, hydrophobicity, conservation, entropy, and pairwise identity. MSATK is designed for command-line use, Python notebooks, teaching, manuscript supplementary outputs, and automated workflows requiring stable filenames and reproducible metadata. The software is freely available under the MIT License and can be installed from source, PyPI, or Conda-compatible environments.

# Introduction

Multiple sequence alignments (MSAs) are a fundamental data structure in bioinformatics. They support phylogenetic inference, molecular evolution, comparative genomics, protein-family analysis, metagenomics, and viral genomic surveillance. Before these downstream analyses, researchers commonly need to evaluate alignment completeness, gap structure, sequence quality, conservation, compositional bias, pairwise similarity, and, for coding sequences, codon-level properties. These summaries are frequently produced using project-specific scripts or by combining several specialized tools, which can make quality-control decisions difficult to reproduce across projects.

Existing software provides excellent support for related but distinct stages of alignment analysis. MAFFT [@katoh2013mafft], MUSCLE [@edgar2004muscle], and Clustal Omega [@sievers2011clustal] construct alignments; Jalview [@waterhouse2009jalview] and AliView [@larsson2014aliview] support interactive alignment inspection; trimAl [@capella2009trimal] assists with alignment trimming; and IQ-TREE [@minh2020iqtree], RAxML [@stamatakis2014raxml], HyPhy [@pond2020hyphy], and PAML [@yang2007paml] support downstream inference. These tools are not intended to provide a single standardized, report-oriented summary of alignment quality, composition, codon usage, warnings, plots, and reproducible output metadata.

MSATK fills this gap by providing a "FastQC-like" profiling layer for multiple sequence alignments. Given an input alignment, MSATK produces a structured output directory containing summary metadata, quality-control warnings, tabular statistics, optional plots, and shareable reports. The goal is not to replace alignment algorithms, manual visualization tools, trimming tools, phylogenetic software, or selection-analysis packages. Instead, MSATK helps users rapidly assess whether an alignment is complete, coherent, codon-aware when expected, and suitable for subsequent analysis.

# Software Overview

MSATK is both a Python package and a command-line application. The primary command profiles an input alignment and writes a complete output directory:

```bash
msatk profile alignment.fasta --out msatk_results/
```

A bundled demonstration command provides a quick smoke test and example report:

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

Specialized Python interfaces are available for codon and protein workflows through `CodonProfiler` and `ProteinProfiler`. `AlignmentProfiler` remains available as an alias for `MSATK`.

# Supported Inputs and Outputs

MSATK supports common multiple-alignment and alignment-derived formats, including FASTA and aligned FASTA, PHYLIP and relaxed PHYLIP, NEXUS matrix blocks, CLUSTAL, Stockholm, A3M, MAF, SAM, and BAM/CRAM-derived summaries through the optional `pysam` dependency. Supported analysis modes include nucleotide, RNA, protein, codon, translated CDS-like, mixed, and auto-detected inputs.

Each profiling run writes a reproducible output directory:

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

Core tables include `alignment_summary.csv`, `per_sequence_stats.csv`, `per_site_stats.csv`, `pairwise_identity_matrix.csv`, `distance_matrix.csv`, `composition_summary.csv`, `entropy_by_site.csv`, `sequence_embeddings.csv`, and `qc_flags.csv`. Codon-aware runs additionally write `codon_usage.csv`, `rscu.csv`, `codon_position_gc.csv`, `codon_position_entropy.csv`, `stop_codon_report.csv`, and `translated_amino_acid_summary.csv`. Protein runs write amino-acid composition, residue-class, and hydrophobicity summaries.

# Implementation

MSATK is implemented in Python and organized as modular layers for input parsing, validation, molecule-type detection, statistics, codon analysis, protein analysis, plotting, embeddings, report generation, and output writing. The internal workflow is:

```text
Input reader -> validation -> detection -> statistics -> plotting -> reporting -> output writer
```

The core package uses lightweight internal parsers for common MSA formats. Optional dependencies extend the user experience: pandas [@mckinney2010pandas] provides notebook-friendly DataFrame outputs, Matplotlib [@hunter2007matplotlib] enables PNG figures, scikit-learn [@pedregosa2011scikit] enables dimensionality-reduction embeddings, and pysam [@pysam] enables BAM and CRAM parsing through HTSlib-compatible interfaces.

The validation layer supports permissive and strict modes. In permissive mode, MSATK warns and continues where scientifically reasonable, for example by padding unequal sequence lengths for site-level summaries. In strict mode, duplicate identifiers, invalid characters, and unequal sequence lengths raise explicit validation errors. Command-line options expose these behaviors through `--permissive` and `--strict`.

Each run records the MSATK version, Python version, platform, input path, detected format, detected molecule type, command-line parameters, timestamp, validation mode, and output schema version. These metadata make the output directory suitable for workflow automation, manuscript supplements, and later audit.

# Alignment Statistics and Quality Control

MSATK computes alignment-wide, per-sequence, and per-site summaries. Alignment-level metrics include sequence count, alignment length, rectangularity, gap fraction, conserved sites, variable sites, parsimony-informative sites, singleton sites, mean pairwise identity, mean pairwise distance, entropy summaries, and GC content for nucleotide-like alignments.

Per-sequence outputs include raw length, ungapped length, gap count, gap fraction, ambiguous-character count, ambiguous-character fraction, GC content, unique residue count, mean identity to other sequences, and a simple outlier score. These values help identify short, gap-rich, ambiguous, compositionally unusual, or low-similarity sequences.

Per-site outputs include dominant character, dominant-character frequency, gap fraction, missing fraction, Shannon entropy [@shannon1948], conservation score, observed-state count, variable-site status, parsimony-informative status, singleton status, and codon position. These summaries support trimming decisions, downstream masking, and inspection of highly variable or poorly covered regions.

Quality-control warnings are written to the terminal, report, `summary.json`, and `qc_warnings.txt`. Warnings include inconsistent sequence lengths, high-gap sequences, high-gap sites, codon alignment lengths not divisible by three, internal stop codons, and unusually low mean identity to the rest of the alignment.

# Codon-Aware and Protein-Specific Summaries

For codon-aware alignments, MSATK reports codon usage, relative synonymous codon usage (RSCU), codon-position GC content, codon-position entropy, stop codons, and translated amino-acid summaries. RSCU summarizes observed codon use relative to expected use within synonymous codon families [@sharp1986rscu]. Frame-related warnings help users detect alignments that may be unsuitable for downstream selection-analysis tools.

For protein alignments, MSATK reports amino-acid composition, residue-class composition, hydrophobicity summaries using the Kyte-Doolittle scale [@kyte1982hydropathy], site-level entropy, conservation, and pairwise identity. These summaries provide rapid context for protein-family alignments before structural, functional, or evolutionary interpretation.

Visualization and reporting functions produce gap profiles, entropy plots, conservation plots, composition plots, pairwise identity heatmaps, pairwise distance distributions, pairwise identity distributions, and sequence embeddings when the relevant optional dependencies are available. HTML and Markdown reports combine these plots with summary statistics, warnings, and plain-language interpretations.

# Example Applications

To demonstrate MSATK on a coding-sequence alignment, a user can run:

```bash
msatk profile tests/data/codon/valid_codon_alignment.fasta --type codon --out codon_results --force
```

MSATK detects a codon-aware nucleotide alignment, reports sequence count and alignment length, computes missingness, summarizes variable sites and mean pairwise identity, and writes codon-specific tables including codon usage, RSCU, codon-position GC, codon-position entropy, stop-codon reports, and translated amino-acid summaries. A related fixture containing an internal stop codon triggers a plain-language warning identifying the affected sequence. These outputs provide a compact pre-analysis check before running selection-analysis tools such as HyPhy or PAML.

The same workflow can be applied to protein alignments:

```bash
msatk protein tests/data/protein/protein_with_gaps.faa --out protein_results --force
```

In this mode, MSATK writes amino-acid composition, residue-class, hydrophobicity, pairwise identity, entropy, and conservation summaries. For workflow developers, all outputs use deterministic filenames under a single output directory, allowing direct integration with Snakemake, Nextflow, or HPC batch jobs.

# Availability and Installation

MSATK is freely available under the MIT License at `https://github.com/aglucaci/MSATK`. The version described in this manuscript is v0.1.0. Source installation and editable development installation are supported with:

```bash
pip install -e ".[dev,all,docs]"
```

MSATK includes a local Conda recipe and a Bioconda submission template. The intended bioinformatics-native installation path is:

```bash
mamba install -c bioconda -c conda-forge msatk
```

Until a Bioconda package is published, users can create a Conda-managed environment from the provided `environment.yml` and install MSATK from source or PyPI. Documentation, synthetic test datasets, GitHub Actions workflows, and release checklists are included with the repository.

# Testing and Reproducibility

The repository includes unit tests, integration tests, regression tests, synthetic alignment fixtures, documentation checks, packaging files, and continuous integration workflows. Test data cover FASTA, aligned FASTA, A3M, PHYLIP, relaxed PHYLIP, CLUSTAL, Stockholm, NEXUS, MAF, SAM, and optional BAM/CRAM fixtures where HTSlib-compatible tooling is available.

MSATK outputs are designed to be reproducible and workflow-friendly. The output directory contains stable filenames, versioned `summary.json`, `parameters.yaml`, `msatk.log`, CSV tables, reports, and optional figures. Continuous integration runs tests across supported Python versions, checks type annotations with mypy, validates style with Ruff, builds documentation, and exercises Conda packaging.

# Discussion

MSATK contributes a standardized profiling and reporting layer for multiple sequence alignments. It lowers the barrier to alignment quality control by turning an input alignment into a coherent set of statistics, warnings, tables, figures, and reports with one command. This is useful for teaching, exploratory data analysis, manuscript supplementary materials, and production bioinformatics workflows.

MSATK complements, rather than replaces, existing tools:

| Tool | Alignment construction | Visualization | QC report | Codon summaries | Reproducible tables |
| --- | ---: | ---: | ---: | ---: | ---: |
| MAFFT | Yes | No | No | No | No |
| Jalview | No | Yes | Limited | Limited | Limited |
| trimAl | No | No | Partial | No | Partial |
| HyPhy | No | No | No | Yes, downstream | Yes |
| MSATK | No | Report-based | Yes | Yes | Yes |

MSATK is not intended to construct alignments, infer phylogenies, perform formal selection analysis, or replace manual inspection of complex alignments. Its main role is the intermediate quality-control and reporting step between alignment generation and downstream biological inference. Current limitations include the computational cost of pairwise summaries for very large alignments and the need for further validation on large pangenome, metagenomic, and surveillance-scale datasets.

Future development will focus on interactive reports, more scalable pairwise computations, richer outlier detection, expanded workflow wrappers, Galaxy integration, containerized deployments, and stronger support for very large alignments.

# Conclusion

MSATK provides a lightweight, reproducible, and user-friendly framework for profiling multiple sequence alignments. By combining alignment statistics, quality-control warnings, codon-aware summaries, protein summaries, visualizations, and standalone reports, MSATK helps researchers rapidly assess alignment quality and prepare data for evolutionary, comparative, and functional analyses.

# Acknowledgements

MSATK builds on the broader open-source scientific Python and bioinformatics ecosystem. The author thanks contributors to Python, pandas, Matplotlib, scikit-learn, pysam, and the alignment and evolutionary-analysis tools that motivate reproducible alignment quality-control workflows.

# Funding

No external funding was received for this work unless otherwise specified by the author.

# Conflict of Interest

The author declares no competing interests.

# References
