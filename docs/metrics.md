# Metric dictionary

Generated from `msatk.stats.registry`. Every field MSATK writes is defined here.

Character states: **residue** (canonical), **ambiguous** (IUPAC/B/Z/J), **missing** (`?`; `N` for nucleotides; `X` for proteins; padding), **gap** (`-`, `.`), **stop** (`*`).
Undefined values are written as empty CSV cells / JSON `null`, never as zero.

## Alignment-level

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `tool` | Name of the generating tool (always 'MSATK'). |  |  |  | all |
| `number_of_sequences` | Number of sequences (rows), N. |  | >= 1 |  | all |
| `alignment_length` | Aligned width L (number of columns). |  | >= 1 |  | all |
| `total_cells` | N x L matrix cells, including padding. | `N * L` |  |  | all |
| `rectangular_alignment` | True when every input sequence has the same aligned length. |  |  |  | all |
| `sequence_lengths_consistent` | Alias of rectangular_alignment (kept for compatibility). |  |  |  | all |
| `padded_cells` | Cells added to shorter input rows to form a rectangle. Padding is recorded as missing data ('?'), never as alignment gaps. |  | >= 0 |  | all |
| `detected_format` | Input file format used for parsing. |  |  |  | all |
| `detected_file_format` | Alias of detected_format. |  |  |  | all |
| `molecule_type` | Molecule type used for analysis (detected or user-specified). |  |  |  | all |
| `detected_molecule_type` | Alias of molecule_type. |  |  |  | all |
| `length_divisible_by_three` | True when L mod 3 == 0. |  |  |  | all |
| `appears_codon_aware` | Nucleotide alignment with L divisible by 3 and no frame warnings. |  |  |  | dna, rna, codon |
| `stop_codons_exist` | Any in-frame stop codon (TAA/TAG/TGA) observed. |  |  |  | dna, rna, codon |
| `stop_codon_count` | Number of in-frame stop codons over all sequences. |  |  |  | dna, rna, codon |
| `internal_stop_codons_exist` | Any stop codon before the final codon. |  |  |  | dna, rna, codon |
| `frameshift_warnings_exist` | Any frame warning (L or ungapped length not divisible by 3, or partially gapped codons). |  |  |  | dna, rna, codon |
| `frameshift_warning_count` | Number of frame warnings. |  |  |  | dna, rna, codon |
| `translated_cds_mode` | True when analysed as translated CDS. |  |  |  | all |
| `duplicate_id_count` | Number of sequence identifiers that occur more than once. |  |  |  | all |
| `unique_sequence_count` | Number of distinct aligned sequences (exact string identity after case normalisation). |  |  |  | all |
| `duplicate_sequence_count` | Sequences identical to an earlier sequence. |  |  |  | all |
| `residue_cells` | Cells holding a canonical residue. *Only canonical residues counted; gaps, missing and ambiguous cells excluded.* |  |  |  | all |
| `ambiguous_cells` | Cells holding a partially informative ambiguity code. |  |  |  | all |
| `missing_cells` | Cells holding missing data (excludes padding). |  |  |  | all |
| `stop_cells` | Cells holding '*' (protein alignments). |  |  |  | all |
| `invalid_cells` | Cells holding characters outside the molecule's alphabet. |  |  |  | all |
| `total_gaps` | Number of gap cells ('-' or '.'). |  |  |  | all |
| `gap_fraction` | Fraction of all cells that are gaps. | `total_gaps / (N*L)` | [0, 1] |  | all |
| `occupancy` | Fraction of all cells holding a canonical residue. | `residue_cells / (N*L)` | [0, 1] |  | all |
| `ambiguous_fraction` | Fraction of non-padding cells holding ambiguity codes. |  | [0, 1] |  | all |
| `mean_sequence_length` | Mean ungapped sequence length. |  |  | No sequences. | all |
| `median_sequence_length` | Median ungapped sequence length. |  |  |  | all |
| `min_sequence_length` | Minimum ungapped sequence length. |  |  |  | all |
| `max_sequence_length` | Maximum ungapped sequence length. |  |  |  | all |
| `sequence_length_sd` | Population standard deviation of ungapped lengths. |  |  |  | all |
| `constant_sites` | Columns with exactly one observed canonical state. All-gap/all-missing columns are NOT constant; they are counted as unobserved. *Gaps, missing (?, N for nucleotides, X for proteins) and ambiguity codes excluded.* |  |  |  | all |
| `conserved_sites` | Alias of constant_sites (kept for compatibility). |  |  |  | all |
| `variable_sites` | Columns with >= 2 observed canonical states. *Gaps, missing (?, N for nucleotides, X for proteins) and ambiguity codes excluded.* |  |  |  | all |
| `unobserved_sites` | Columns with no canonical residue (statistics undefined). |  |  |  | all |
| `parsimony_informative_sites` | Columns with >= 2 states each observed in >= 2 sequences. *Gaps, missing (?, N for nucleotides, X for proteins) and ambiguity codes excluded.* |  |  |  | all |
| `singleton_sites` | Variable columns that are not parsimony-informative (MEGA definition). *Gaps, missing (?, N for nucleotides, X for proteins) and ambiguity codes excluded.* |  |  |  | all |
| `gap_free_sites` | Columns without any gap. |  |  |  | all |
| `complete_sites` | Columns in which every sequence carries a canonical residue. |  |  |  | all |
| `entropy_mean` | Mean site Shannon entropy (bits) over observed columns. *Unobserved columns excluded.* |  |  | No observed column. | all |
| `entropy_median` | Median site entropy over observed columns. |  |  |  | all |
| `normalized_entropy_mean` | Mean of entropy / log2(K) over observed columns. |  | [0, 1] |  | all |
| `mean_consensus_identity` | Mean over sequences of identity to the majority consensus. |  | [0, 1] |  | all |
| `pairwise_method` | How pairwise identity was computed: exact (all pairs, full matrix), all_pairs (all pairs, no matrix), sampled (deterministic random sample with replacement) or none. |  |  |  | all |
| `pairwise_pairs_total` | N(N-1)/2. |  |  |  | all |
| `pairwise_pairs_evaluated` | Number of pairs actually compared. |  |  |  | all |
| `pairwise_pairs_without_comparable_sites` | Evaluated pairs sharing no comparable site (identity undefined, excluded from means). |  |  |  | all |
| `pairwise_sample_seed` | Random seed used for sampled pairwise comparison. |  |  |  | all |
| `mean_pairwise_identity` | Mean pairwise identity over evaluated pairs with >= 1 comparable site. *A site is comparable only when both sequences carry a canonical residue.* | `matches / comparable_sites per pair` | [0, 1] | Pairwise comparison disabled, N < 2, or no pair has comparable sites. | all |
| `median_pairwise_identity` | Median pairwise identity over evaluated pairs. |  |  |  | all |
| `min_pairwise_identity` | Minimum pairwise identity over evaluated pairs. |  |  |  | all |
| `max_pairwise_identity` | Maximum pairwise identity over evaluated pairs. |  |  |  | all |
| `mean_pairwise_distance` | 1 - mean_pairwise_identity (uncorrected p-distance). |  |  |  | all |
| `gc_content` | Pooled G+C fraction over all canonical bases. *Only canonical residues counted; gaps, missing and ambiguous cells excluded.* | `(G + C) / (A + C + G + T)` |  | No canonical bases (always null for proteins). | dna, rna, codon |
| `gc_skew` | (G - C) / (G + C) over all canonical bases. |  | [-1, 1] |  | dna, rna, codon |
| `at_skew` | (A - T) / (A + T) over all canonical bases. |  | [-1, 1] |  | dna, rna, codon |
| `popgen_sites_used` | Number of complete sites used for population-genetic statistics. |  |  |  | dna, rna, codon |
| `segregating_sites` | Number of polymorphic complete sites, S. *Computed on complete sites (every sequence carries a canonical base). Assumes homologous population samples; interpret cautiously for divergent multi-species alignments.* |  |  |  | dna, rna, codon |
| `nucleotide_diversity_pi` | Nucleotide diversity: mean pairwise differences per complete site (Nei & Li 1979). *Computed on complete sites (every sequence carries a canonical base). Assumes homologous population samples; interpret cautiously for divergent multi-species alignments.* | `sum_sites [C(N,2) - sum_c C(n_c,2)] / C(N,2) / sites_used` |  | N < 2 or no complete site. | dna, rna, codon |
| `theta_watterson` | Watterson's theta per site. *Computed on complete sites (every sequence carries a canonical base). Assumes homologous population samples; interpret cautiously for divergent multi-species alignments.* | `S / a1 / sites_used, a1 = sum_{i=1}^{N-1} 1/i` |  |  | dna, rna, codon |
| `tajimas_d` | Tajima's D (Tajima 1989). *Computed on complete sites (every sequence carries a canonical base). Assumes homologous population samples; interpret cautiously for divergent multi-species alignments.* |  |  | N < 4 or S = 0. | dna, rna, codon |
| `haplotype_count` | Distinct haplotypes over complete sites. *Computed on complete sites (every sequence carries a canonical base). Assumes homologous population samples; interpret cautiously for divergent multi-species alignments.* |  |  |  | dna, rna, codon |
| `haplotype_diversity` | Nei's haplotype diversity over complete sites. *Computed on complete sites (every sequence carries a canonical base). Assumes homologous population samples; interpret cautiously for divergent multi-species alignments.* | `N/(N-1) * (1 - sum_i p_i^2)` | [0, 1] |  | dna, rna, codon |
| `transitions_pairwise` | Transition differences summed over all pairs at complete sites. |  |  |  | dna, rna, codon |
| `transversions_pairwise` | Transversion differences summed over all pairs at complete sites. |  |  |  | dna, rna, codon |
| `ts_tv_ratio` | Ratio of pooled pairwise transition to transversion differences (uncorrected). |  |  | No transversions observed. | dna, rna, codon |
| `msatk_version` | MSATK version that produced the output. |  |  |  | all |
| `output_schema_version` | Version of the MSATK output schema. |  |  |  | all |
| `input` | Input path as given. |  |  |  | all |
| `input_sha256` | SHA-256 digest of the input file. |  |  |  | all |
| `generated_at` | UTC timestamp of the run (ISO 8601). |  |  |  | all |
| `python_version` | Python interpreter version. |  |  |  | all |
| `platform` | Operating system / platform string. |  |  |  | all |
| `coordinate_convention` | Coordinate convention of site indices (1-based, inclusive). |  |  |  | all |
| `interpretation` | Plain-language interpretation of the profile. |  |  |  | all |
| `recommended_next_steps` | Suggested follow-up actions. |  |  |  | all |

## Sequence-level

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `row_index` | 1-based row index in the input alignment. |  |  |  | all |
| `sequence_id` | Sequence identifier. |  |  |  | all |
| `original_length` | Length of the row as provided (including gaps), before padding. |  |  |  | all |
| `aligned_width` | Alignment width L. |  |  |  | all |
| `padded_cells` | Missing-data cells added to reach width L. |  |  |  | all |
| `ungapped_length` | Non-gap, non-padding characters in the row. |  |  |  | all |
| `residue_count` | Canonical residues in the row. |  |  |  | all |
| `span_start` | First non-gap column (1-based). |  |  | All-gap row. | all |
| `span_end` | Last non-gap column (1-based). |  |  | All-gap row. | all |
| `span_length` | span_end - span_start + 1 (0 for an all-gap row). |  |  |  | all |
| `gap_count` | Gap cells in the row. |  |  |  | all |
| `gap_fraction` | gap_count / L. |  | [0, 1] |  | all |
| `leading_gaps` | Gap cells before the first non-gap character. |  |  |  | all |
| `trailing_gaps` | Gap cells after the last non-gap character. |  |  |  | all |
| `internal_gap_count` | Gap cells between span_start and span_end. |  |  |  | all |
| `gap_run_count` | Number of contiguous gap runs (including terminal runs). |  |  |  | all |
| `max_gap_run` | Length of the longest gap run. |  |  |  | all |
| `missing_count` | Missing-data cells in the row (padding excluded). |  |  |  | all |
| `ambiguous_count` | Ambiguity-code cells in the row. |  |  |  | all |
| `ambiguous_fraction` | ambiguous_count / L. |  |  |  | all |
| `occupancy` | residue_count / L. |  | [0, 1] |  | all |
| `unique_residue_count` | Distinct canonical residues observed in the row. |  |  |  | all |
| `consensus_comparable_sites` | Columns where the row has a residue and the consensus is defined. |  |  |  | all |
| `consensus_identity` | Fraction of comparable columns where the row matches the majority consensus. |  | [0, 1] | No comparable column. | all |
| `composition_distance` | Total-variation distance between the row's residue composition and the pooled composition. | `0.5 * sum_c |p_c(row) - p_c(alignment)|` | [0, 1] |  | all |
| `duplicate_of` | Identifier of the first earlier row with an identical sequence, if any. |  |  |  | all |
| `gc_content` | G+C fraction among canonical bases. |  |  |  | dna, rna, codon |
| `gc_skew` | (G - C) / (G + C). |  |  |  | dna, rna, codon |
| `at_skew` | (A - T) / (A + T). |  |  |  | dna, rna, codon |
| `mean_identity_to_others` | Mean identity to every other sequence (exact pairwise mode only). |  |  | Pairwise matrix not computed or no comparable partner. | all |
| `nearest_neighbor` | Most similar other sequence (exact pairwise mode only). |  |  |  | all |
| `nearest_neighbor_identity` | Identity to the nearest neighbour. |  |  |  | all |
| `outlier_score` | 1 - consensus_identity. |  | [0, 1] |  | all |
| `consensus_identity_robust_z` | Robust z-score of consensus_identity: (x - median) / (1.4826 * MAD). |  |  | Fewer than 3 defined values. | all |
| `composition_robust_z` | Robust z-score of composition_distance. |  |  |  | all |

## Site-level

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `site_index` | 1-based column index in the input alignment. |  |  |  | all |
| `site_class` | constant, variable or unobserved (no canonical residue). |  |  |  | all |
| `observed_count` | Sequences with a canonical residue in the column. *Only canonical residues counted; gaps, missing and ambiguous cells excluded.* |  |  |  | all |
| `occupancy` | observed_count / N. |  | [0, 1] |  | all |
| `gap_count` | Gap cells in the column. |  |  |  | all |
| `gap_fraction` | gap_count / N. |  | [0, 1] |  | all |
| `missing_fraction` | Missing-data cells / N (gaps and ambiguity codes not included). |  | [0, 1] |  | all |
| `ambiguous_fraction` | Ambiguity-code cells / N. |  | [0, 1] |  | all |
| `dominant_character` | Most frequent canonical residue (ties broken by alphabet order). |  |  | Unobserved column. | all |
| `dominant_character_frequency` | Frequency of the dominant residue among observed residues. | `max_c n_c / observed_count` | [1/K, 1] | Unobserved column. | all |
| `conservation_score` | Alias of dominant_character_frequency. |  |  | Unobserved column. | all |
| `consensus_support` | Fraction of ALL sequences carrying the dominant residue (penalises gaps). | `max_c n_c / N` | [0, 1] |  | all |
| `observed_state_count` | Number of distinct canonical residues in the column. |  |  |  | all |
| `entropy` | Shannon entropy in bits of the canonical residue distribution. *Gaps, missing (?, N for nucleotides, X for proteins) and ambiguity codes excluded.* | `-sum_c p_c log2 p_c` | [0, log2 K] | Unobserved column. | all |
| `normalized_entropy` | entropy / log2(K), K = canonical alphabet size. |  | [0, 1] |  | all |
| `effective_states` | Effective number of states, 2^entropy. |  | [1, K] |  | all |
| `information_content` | log2(K) - entropy (no small-sample correction). |  | [0, log2 K] |  | all |
| `constant` | True when exactly one canonical state is observed. |  |  |  | all |
| `variable` | True when >= 2 canonical states are observed. |  |  |  | all |
| `parsimony_informative` | True when >= 2 states each occur in >= 2 sequences. |  |  |  | all |
| `singleton` | True when variable but not parsimony-informative. |  |  |  | all |
| `gc_content` | G+C fraction among canonical bases in the column. |  |  |  | dna, rna, codon |
| `codon_position` | Codon position (1, 2, 3) of the column. |  |  |  | codon |

## Pair-level

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_a` | First sequence of the pair. |  |  |  | all |
| `sequence_b` | Second sequence of the pair. |  |  |  | all |
| `comparable_sites` | Columns where both sequences carry a canonical residue. |  |  |  | all |
| `mismatches` | Comparable columns with different residues. |  |  |  | all |
| `identity` | (comparable_sites - mismatches) / comparable_sites. |  |  | No comparable sites. | all |
| `p_distance` | 1 - identity (uncorrected). |  |  | No comparable sites. | all |

