# Outputs

Default profile output:

```text
alignment_msatk/
├── report.html
├── report.md
├── summary.json
├── parameters.yaml
├── msatk.log
├── qc_warnings.txt
├── tables/
│   ├── alignment_summary.csv
│   ├── per_sequence_stats.csv
│   ├── per_site_stats.csv
│   ├── pairwise_identity_matrix.csv
│   ├── distance_matrix.csv
│   ├── composition_summary.csv
│   ├── entropy_by_site.csv
│   ├── sequence_embeddings.csv
│   └── qc_flags.csv
└── figures/
```

Codon alignments additionally include codon usage, RSCU, codon-position GC, codon-position entropy, translated amino-acid summaries, and stop codon reports.

## Stable Table Schemas

`per_sequence_stats.csv` columns:

- `sequence_id`
- `raw_length`
- `ungapped_length`
- `gap_count`
- `gap_fraction`
- `ambiguous_count`
- `ambiguous_fraction`
- `gc_content`
- `unique_residue_count`
- `mean_identity_to_others`
- `outlier_score`

`per_site_stats.csv` columns:

- `site_index`
- `dominant_character`
- `dominant_character_frequency`
- `gap_fraction`
- `missing_fraction`
- `entropy`
- `conservation_score`
- `observed_state_count`
- `variable`
- `parsimony_informative`
- `singleton`
- `codon_position`

Every run writes `summary.json`, `parameters.yaml`, and `msatk.log` with MSATK version, Python version, platform, detected format, detected molecule type, run parameters, and output schema version.
