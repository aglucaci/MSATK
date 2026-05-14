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

## Supported Input Formats

MSATK supports common multiple-sequence and read-alignment formats:

- FASTA and aligned FASTA
- A3M
- PHYLIP and relaxed PHYLIP
- CLUSTAL
- Stockholm
- NEXUS matrix blocks
- MAF
- SAM
- BAM
- CRAM

SAM/BAM/CRAM are treated as alignment-derived inputs: MSATK uses mapped reads and CIGAR operations to build profileable aligned records. SAM is parsed directly. BAM and CRAM require the optional `pysam` dependency.

## Auto-Detection Fields

MSATK automatically infers and writes these fields to `summary.json`, `tables/alignment_summary.csv`, and the report:

- `detected_file_format`
- `detected_molecule_type`
- `alignment_length`
- `sequence_lengths_consistent`
- `length_divisible_by_three`
- `appears_codon_aware`
- `stop_codons_exist`
- `stop_codon_count`
- `internal_stop_codons_exist`
- `frameshift_warnings_exist`
- `frameshift_warning_count`
- `translated_cds_mode`
