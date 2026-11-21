# Outputs

`msatk profile` writes a predictable directory (output schema version **2.0**):

```text
alignment_msatk/
├── report.html              # interactive offline dashboard
├── report.md
├── summary.json             # alignment-level statistics + provenance
├── metric_dictionary.json   # definition of every field
├── qc_findings.json
├── qc_warnings.txt
├── parameters.yaml          # complete run configuration
├── run_manifest.json        # configuration, dependency versions, list of outputs
├── msatk.log
├── assets/
├── tables/
│   ├── alignment_summary.csv
│   ├── per_sequence_stats.csv
│   ├── per_site_stats.csv
│   ├── composition_summary.csv
│   ├── entropy_by_site.csv
│   ├── gap_blocks.csv
│   ├── qc_findings.csv
│   ├── qc_flags.csv
│   ├── pairwise_identity_matrix.csv   # exact pairwise mode only
│   ├── distance_matrix.csv            # exact pairwise mode only
│   ├── pairwise_comparisons.csv       # exact pairwise mode only
│   └── sequence_embeddings.csv        # only with --embeddings
├── figures/                 # one file per requested format (png/svg/pdf/tiff)
│   ├── alignment_overview.*
│   ├── sequence_logo.*
│   └── ...
└── ml/                      # only with --encodings
```

Specialized analyses add tables when they apply (see
[Alignment types](alignment_types.md)): `site_substitution_scores.csv` (proteins),
`rna_base_pairs.csv` (Stockholm `SS_cons`), `profile_states.csv` (A2M/A3M) and
`partitions.csv` (charsets or `--partitions`). Their summaries are nested under
`summary.json → specialized`, and `analysis_applicability` records why each analysis ran or
was skipped.

Block, read, mapping, graph and profile-HMM inputs produce a smaller directory:
`summary.json` (statistics and provenance), `tables/*.csv` (for example `blocks.csv` and
`genomes.csv`; `references.csv` and `mapq_histogram.csv`; `queries.csv` and `targets.csv`;
`segments.csv` and `paths.csv`; `match_states.csv`), `metric_dictionary.json`,
`report.md`, `report.html` and `run_manifest.json`. Intervals in these tables are 0-based
and half-open.

Codon alignments add codon usage, RSCU, codon-position GC and entropy, translated
amino-acid summaries and stop-codon reports. Protein alignments add amino-acid,
residue-class and hydrophobicity summaries.

## Conventions

* Site coordinates are **1-based, inclusive** and refer to the input alignment.
* Undefined values (e.g. entropy of an all-gap column, identity of non-overlapping
  sequences) are empty CSV cells / JSON `null` — never zero.
* Every column is defined in the [metric dictionary](metrics.md); the same
  definitions ship with each run as `metric_dictionary.json`.

## Provenance

`summary.json` and `run_manifest.json` record the input SHA-256, MSATK version,
timestamp, Python version, platform, all parameters (including the pairwise policy
and random seed) and installed dependency versions.

## Changes from schema 1.0

* `raw_length` was removed (it was neither the original nor the ungapped length).
  Use `original_length`, `ungapped_length`, `aligned_width` and `span_*`.
* `conserved_sites` no longer counts all-gap columns; see `unobserved_sites`.
* `mean_identity_to_others` is null unless an exact pairwise matrix is computed;
  `outlier_score` is now `1 - consensus_identity` (O(N·L), always available).
* `codon_position` appears only for codon alignments; `gc_content` only for nucleotides.
* Pairwise matrices and embeddings are no longer written unconditionally.
