# Changelog

## 0.2.0 - 2026-11-20

Scientific semantics (breaking output-schema changes; schema version 2.0):

- New molecule-aware character policy (`msatk.core.alphabet`) shared by all metrics.
  Fixes asparagine (`N`) and `X` being treated as missing data in protein alignments.
- All-gap/all-missing columns are now `unobserved`, not constant/conserved; their
  entropy and conservation are `null`.
- Pairwise identity is `null` (not 0) when two sequences share no comparable site;
  such pairs are counted and excluded from averages.
- Unequal-length input is padded with missing data (`?`), never with gaps, and the
  number of padded cells is reported.
- Singleton sites follow the standard definition (variable, not parsimony-informative).
- `codon_position` is only reported for codon alignments; `gc_content` only for nucleotides.
- Sequence lengths are split into `original_length`, `aligned_width`, `ungapped_length`,
  `residue_count` and `span_*`; `raw_length` was removed.

Statistics:

- New vectorised `msatk.stats` package (alignment, sequence, site, pairwise,
  information, diversity) and a metric dictionary (`msatk.stats.registry`,
  `metric_dictionary.json`, `docs/metrics.md`) defining every field.
- New fields: occupancy, ambiguous/missing/stop/invalid cells, duplicate IDs and
  sequences, complete/gap-free/unobserved sites, normalised entropy, effective
  states, information content, consensus support, gap runs and spans, GC/AT skew,
  consensus identity, composition distance, robust z-scores, nearest neighbour,
  gap blocks, and population-genetic statistics (π, Watterson's θ, Tajima's D,
  haplotype diversity, Ts/Tv) validated against scikit-allel.
- Pairwise identity uses blocked matrix products with a size-aware policy
  (`--pairwise auto|exact|sampled|none`).
- Frame/stop-codon detection vectorised (~7x faster on large alignments).

QC:

- Structured, severity-ranked findings with codes, affected items and suggestions
  (`qc_findings.json/.csv`); new rules for duplicates, invalid characters,
  unobserved sites, composition outliers, ambiguity, non-overlapping pairs,
  partially gapped codons and internal protein stops.

Visualization and reports:

- Residue-level alignment figures with consensus and annotation tracks, row ordering,
  difference highlighting and custom themes (`msatk plot`); sequence logos (`msatk logo`).
- Figures are written in every requested format (PNG/SVG/PDF/TIFF); unsupported
  formats are rejected.
- New interactive offline HTML report: residue viewer, site tracks, sortable/searchable
  tables with CSV export, pairwise heatmap, embedding plot, provenance and metric
  dictionary; `--self-contained` single-file mode.

Machine learning:

- New `msatk.encoding` package: lossless categorical, one-hot (missing/fractional/
  separate ambiguity policies), physicochemical and feature-channel tensors;
  coordinate-preserving crops and windows with padding masks; RGB and label images;
  NPZ/PNG exports with `manifest.json` (`msatk encode`, `msatk image`).
- `msatk.ml.AlignmentTensorDataset` plus PyTorch/TensorFlow adapters.
- Embeddings are opt-in and never silently fall back to another method; PCA no
  longer requires scikit-learn.

CLI/API and packaging:

- New commands: `stats`, `plot`, `logo`, `encode`, `image`, `metrics`; new profile
  options `--preset`, `--pairwise`, `--embeddings`, `--encodings`, `--self-contained`.
- New API methods: `sites()`, `sequences()`, `pairwise()`, `qc()`, `gap_blocks()`,
  `plot_alignment()`, `plot_logo()`, `to_tensor()`, `to_image()`, `export_encoding()`,
  `report()`. Earlier method names remain.
- Provenance: input SHA-256, dependency versions and full configuration in
  `summary.json`, `parameters.yaml` and `run_manifest.json`.
- NumPy is now a core dependency; unused `seaborn`, `plotly`, `rich` and `typer`
  extras were removed.
- Benchmark harness (`benchmarks/run_benchmarks.py`).

## 0.1.0 - 2026-05-14

- Initial MSATK package and CLI.
- Added `msatk profile`, `qc`, `codon`, `protein`, `embed`, and `report`.
- Added `msatk demo`.
- Added FASTA/A3M, PHYLIP, CLUSTAL, Stockholm, and NEXUS-style readers.
- Added alignment summaries, per-sequence stats, per-site stats, pairwise identity, distance matrices, QC warnings, composition tables, codon usage, RSCU, protein residue-class summaries, embeddings, plots, and HTML/Markdown reports.
- Added FastQC-style root outputs: `report.html`, `summary.json`, `qc_warnings.txt`, `parameters.yaml`, and `msatk.log`.
- Added plain-language interpretation and recommended next steps.
