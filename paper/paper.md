---
title: "MSATK: comprehensive statistical profiling, visualization and machine-learning representations of multiple sequence alignments"
tags:
  - Python
  - bioinformatics
  - multiple sequence alignment
  - quality control
  - phylogenetics
  - machine learning
authors:
  - name: Alexander Lucaci
    affiliation: 1
affiliations:
  - name: MSATK Project
    index: 1
date: 20 November 2026
bibliography: references.bib
---

<!--
DRAFT — Bioinformatics Application Note (target ≤ 2,000 words + 1 figure).
Sections marked [TODO] require experiments that have not yet been run; do not
submit until they are filled with measured results. Word budget per section is
given in comments.
-->

# Abstract

<!-- ~150 words -->

**Summary:** Multiple sequence alignments (MSAs) underpin phylogenetics, molecular
evolution and protein-family analysis, yet the checks applied to them before
downstream inference are usually ad hoc scripts whose treatment of gaps, missing
data and ambiguity is undocumented. MSATK is a Python toolkit and command-line
application that turns an existing alignment into explicitly defined statistics at
alignment, sequence, site and pair resolution, severity-ranked quality-control
findings, publication-quality figures, an interactive offline HTML report, and
reproducible numerical and image representations for machine learning. All metrics
share one molecule-aware character policy, undefined quantities are reported as
missing rather than zero, and every output field is defined in a versioned metric
dictionary. A vectorised core profiles a 5,000 × 10,000 alignment in about five
seconds.

**Availability and implementation:** MSATK is implemented in Python (≥ 3.10) and
released under the MIT licence at https://github.com/aglucaci/MSATK, with
documentation, tests and benchmarks. Version 0.2.0 is described.

# 1 Introduction

<!-- ~300 words -->

Before an MSA is used for tree inference, selection analysis or model training,
analysts need to know how complete it is, where gaps and ambiguity concentrate,
how conserved and informative its columns are, and which sequences are outliers.
Existing tools address parts of this problem: AMAS computes summary statistics and
manipulates alignments [@borowiec2016amas], AliStat reports completeness
[@wong2020alistat], trimAl filters columns [@capella2009trimal], CIAlign cleans and
visualises alignments [@tumescheit2022cialign], and Jalview and AliView support
interactive inspection [@waterhouse2009jalview; @larsson2014aliview].

Three gaps remain. First, conventions are implicit: whether an all-gap column
counts as conserved, whether non-overlapping sequences have identity zero or no
identity, and whether `N` is missing data (DNA) or asparagine (protein) differ
between scripts and change results. Second, statistics, figures and reports are
produced by separate tools with separate parsers. Third, ML studies that encode
alignments as tensors re-implement encodings without recording vocabulary, channel
semantics or coordinate provenance.

MSATK addresses these with a single analysis engine whose conventions are explicit,
tested and exported alongside every result.

# 2 Implementation

<!-- ~450 words -->

MSATK reads FASTA/A3M, PHYLIP, CLUSTAL, Stockholm, NEXUS, MAF and SAM (BAM/CRAM via
pysam [@pysam]) and detects molecule type. The alignment is vectorised once into
three N × L arrays — normalised characters, cell state and canonical-residue index —
under a molecule-specific policy that assigns each cell to *residue*, *ambiguous*,
*missing*, *gap* or *stop*. Unequal-length input is padded with missing data, never
with gaps, and the padding is reported.

All statistics are NumPy [@harris2020numpy] array operations on this shared
encoding. Column-level counts drive entropy [@shannon1948], information content,
conservation, parsimony-informative and singleton classification; columns without
residues are classified as *unobserved* and receive null statistics. Sequence-level
statistics include gap-run structure, spans, ambiguity, composition, identity to
the majority consensus and composition distance, with robust (median/MAD) z-scores
for outlier detection that avoid the O(N²) cost of pairwise comparison.

Pairwise identity is computed over mutually comparable sites as Σ_c X_c X_cᵀ over
one-hot indicator matrices in row blocks. A size-aware policy chooses exact
comparison up to 1,000 sequences and a seeded sample of pairs above that; the method
and seed are recorded. For nucleotide alignments, MSATK computes segregating sites,
nucleotide diversity [@nei1979], Watterson's θ [@watterson1975], Tajima's D
[@tajima1989], haplotype diversity and pooled Ts/Tv on complete sites.

Each statistic is registered with its resolution, formula, denominator,
missing-data policy, range, applicability and the condition under which it is
undefined; the registry is written as `metric_dictionary.json` with every run.
Provenance (input SHA-256, version, parameters, seed, dependency versions,
coordinate convention) is written to `summary.json`, `parameters.yaml` and
`run_manifest.json`.

The CLI (`msatk profile|stats|qc|plot|logo|encode|image|...`) and the Python API
(`AlignmentProfiler.summary()`, `.sites()`, `.sequences()`, `.pairwise()`, `.qc()`,
`.plot_alignment()`, `.to_tensor()`, `.report()`) call the same cached engine.

# 3 Functionality

<!-- ~450 words -->

**Quality control.** Sixteen explicit rules produce findings with a code,
severity, explanation, suggested action and the affected sequences or sites (e.g.
high-gap rows or columns, divergent or compositionally atypical sequences,
non-overlapping pairs, internal stop codons, partially gapped codons).

**Codon and protein modules.** Codon alignments add codon usage, RSCU
[@sharp1986rscu], GC1/GC2/GC3, codon-position entropy and stop-codon reports;
protein alignments add amino-acid and residue-class composition and Kyte–Doolittle
hydropathy [@kyte1982hydropathy].

**Visualization.** Residue-level alignment figures show a region with consensus row,
conservation/entropy/gap tracks, row ordering and difference highlighting; sequence
logos use the small-sample-corrected information content [@schneider1986]. Themes
(including a colour-vision-deficiency-safe palette) are configurable via JSON/YAML,
and figures are written as PNG, SVG, PDF or TIFF via Matplotlib [@hunter2007matplotlib].

**Interactive report.** `report.html` requires no server or network access. It
embeds the data as JSON and provides a canvas residue viewer (scroll, zoom, ID
search, jump-to-site, colouring by residue or difference from consensus), clickable
site tracks, sortable and filterable tables with CSV export, a pairwise heatmap,
QC findings, provenance and the metric dictionary; a single-file mode embeds figures.

**ML representations.** MSATK keeps human-facing RGB images separate from faithful
encodings: a categorical encoding over a full IUPAC/amino-acid vocabulary that is
exactly reversible; one-hot channels with explicit gap and missing channels and a
choice of ambiguity policy (missing, fractional IUPAC expansion, or separate
channel); physicochemical channels for proteins [@kyte1982hydropathy;
@zamyatnin1972]; and per-cell QC feature channels. Instead of resizing, alignments
are cut into integer-exact windows with padding masks and per-tile original
coordinates. Every export carries a manifest, and a framework-agnostic dataset class
feeds PyTorch or TensorFlow directly.

# 4 Results

<!-- ~450 words -->

**Correctness.** The test suite (93 tests) checks hand-calculated values for
edge cases including all-gap columns, non-overlapping pairs, padding, molecule-specific
interpretation of `N`, entropy bounds, parsimony and singleton definitions and gap
runs. Nucleotide diversity, Watterson's θ and Tajima's D agree with scikit-allel to
floating-point precision on randomly generated populations. Categorical encodings
reconstruct the input alignment exactly, and one-hot cells sum to one under every
ambiguity policy.

**Performance.** On an Apple M1 Pro (16 GB), computing all statistics without
pairwise comparison took 0.23 s for 100 × 10,000, 1.0 s for 1,000 × 10,000 and
4.4 s (1.2 GB peak memory) for 5,000 × 10,000 alignments, scaling linearly with the
number of cells. Exact all-vs-all identity for 5,000 sequences (12.5 million pairs)
added about 4 s; the default sampled policy added 0.3 s (`docs/benchmarks.md`).

**Comparison with existing tools.** [TODO: agreement of shared metrics — alignment
length, missing-data proportion, variable and parsimony-informative sites — with
AMAS and AliStat on public DNA and protein alignments; runtime on the same inputs;
feature-coverage table versus AMAS, AliStat, CIAlign, trimAl and Jalview, moved to
the Supplement.]

**Use cases.** [TODO: (i) a viral coding-sequence alignment: frame integrity, codon
usage, positional entropy, outlier sequences and the QC report before selection
analysis; (ii) a protein-family alignment: physicochemical conservation, clustering,
alignment figures and windowed tensor export loaded by PyTorch.]

**Figure 1.** [TODO: four panels built from real MSATK outputs — (A) architecture;
(B) statistical outputs/site tracks; (C) residue-level alignment figure with logo;
(D) ML encodings and windows with manifest.]

# 5 Discussion

<!-- ~150 words -->

MSATK's contribution is not any single statistic but a reproducible
alignment-to-analysis layer in which conventions are explicit, outputs are defined
and machine-readable, and the same engine drives statistics, figures, reports and
ML encodings. MSATK does not build, trim or edit alignments, and it does not perform
model-based evolutionary inference; it complements aligners, trimmers and
phylogenetic software. Population-genetic statistics assume homologous population
samples and use complete sites only. The interactive viewer is limited to alignments
of up to three million cells, and alignments much larger than memory require
streaming, which is planned together with coevolution statistics and alternative
genetic codes.

# Availability and implementation

MSATK v0.2.0 is available at https://github.com/aglucaci/MSATK under the MIT licence
and installable with `pip install msatk`; documentation, the metric dictionary,
tests and the benchmark harness are included in the repository.

# Funding

No external funding was received for this work unless otherwise specified by the author.

# Conflict of Interest

The author declares no competing interests.

# References
