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
| `genetic_code` | NCBI translation table used for stop-codon and frame checks. |  | NCBI table id (1-33) | Protein or non-nucleotide alignment. | dna, rna, codon |
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
| `alignment_kind` | Data model of the input (sequence, block, read, mapping, graph, profile_model). |  |  |  | all |
| `alignment_scope` | single, pairwise (N = 2) or multiple (N ≥ 3). |  |  |  | all |
| `alignment_method` | Alignment method stated by the input format (global/local), else 'unknown'; never inferred. |  |  |  | all |
| `annotations_found` | Format annotations preserved from the input (e.g. GC:SS_cons, insert_states, partitions). |  |  |  | all |
| `specialized` | Nested summaries of the specialized analyses that apply to this input (see the levels below). |  |  |  | all |
| `analysis_applicability` | For each specialized analysis, whether it ran and why or why not. |  |  |  | all |
| `characterization` | Nested summaries of the characterization analyses that ran (levels below). |  |  |  | all |
| `characterization_profile` | Metrics profile used (core, comprehensive, evolutionary, coevolution, all, none or custom). |  |  |  | all |
| `characterization_applicability` | Why each characterization analysis ran or was skipped. |  |  |  | all |

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
| `composition_distance` | Total-variation distance between the row's residue composition and the pooled composition. | `0.5 * sum_c \|p_c(row) - p_c(alignment)\|` | [0, 1] |  | all |
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

## Pairwise alignment (N = 2)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `alignment_method` | Method stated by the input (EMBOSS needle/stretcher = global, water/matcher = local). |  |  |  | all |
| `alignment_columns` | Alignment length in columns. |  |  |  | all |
| `aligned_pairs` | Columns where neither sequence has a gap. |  |  |  | all |
| `identical_pairs` | Aligned pairs of identical canonical residues. |  |  |  | all |
| `mismatched_pairs` | Aligned pairs of different canonical residues. |  |  |  | all |
| `uncertain_pairs` | Aligned pairs involving an ambiguity or missing-data symbol (not scored). |  |  |  | all |
| `identity` | Identical / (identical + mismatched) pairs. |  |  | No comparable pair. | all |
| `identity_over_columns` | Identical pairs / alignment columns (EMBOSS 'Identity'). |  |  |  | all |
| `similarity` | (Identical + BLOSUM62-positive mismatches) / compared pairs. |  |  | Nucleotides. | protein, translated_cds |
| `length_1` | Residues in sequence 1. |  |  |  | all |
| `length_2` | Residues in sequence 2. |  |  |  | all |
| `coverage_1` | Fraction of sequence-1 residues aligned to a residue of sequence 2. |  |  |  | all |
| `coverage_2` | Fraction of sequence-2 residues aligned to a residue of sequence 1. |  |  |  | all |
| `gap_openings` | Internal gap runs in either sequence (terminal gaps excluded). |  |  |  | all |
| `gap_extensions` | Internal gap positions beyond the first of each run. |  |  |  | all |
| `internal_gap_positions` | Gap positions not touching either end. |  |  |  | all |
| `terminal_gap_positions` | Gap positions in leading/trailing runs (end gaps). |  |  |  | all |
| `score_recomputed` | Score under scoring_scheme (end gaps free, EMBOSS default). |  |  |  | all |
| `scoring_scheme` | Scoring used for score_recomputed. |  |  |  | all |
| `score_reported` | Score stated in the input file (EMBOSS), for comparison. |  |  |  | all |
| `program_reported` | Aligner named in the input file. |  |  |  | all |

## Protein substitution scores (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `mean_column_blosum62` | Mean over columns of the mean BLOSUM62 score of all residue pairs in the column. |  |  |  | protein, translated_cds |
| `median_column_blosum62` | Median of the column mean BLOSUM62 scores. |  |  |  | protein, translated_cds |
| `columns_scored_blosum62` | Columns with at least two canonical residues. |  |  |  | protein, translated_cds |

## Protein substitution scores (per site)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `site_index` | 1-based alignment column. |  |  |  | all |
| `residue_pairs` | Unordered pairs of canonical residues in the column. |  |  |  | all |
| `mean_blosum62` | Mean BLOSUM62 score over those pairs. | `(cᵀSc − Σ c_a S_aa) / (n(n−1))` |  | Fewer than two residues. | all |

## RNA secondary structure (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `consensus_base_pairs` | Base pairs in #=GC SS_cons (WUSS / dot-bracket, pseudoknot letters). |  |  |  | dna, rna, codon |
| `consensus_paired_sites` | Columns paired in the consensus structure. |  |  |  | dna, rna, codon |
| `mean_canonical_pair_fraction` | Mean over pairs of the fraction of complete sequences forming a Watson–Crick or GU pair. |  |  |  | dna, rna, codon |
| `base_pairs_fully_canonical` | Pairs canonical in every complete sequence. |  |  |  | dna, rna, codon |
| `base_pairs_with_compensatory_change` | Pairs where at least one sequence has a double (compensatory) substitution that preserves pairing. |  |  |  | dna, rna, codon |
| `compensatory_pair_observations` | Sequence × pair observations that are compensatory. |  |  |  | dna, rna, codon |
| `consistent_pair_observations` | Sequence × pair observations that are consistent (single change preserving pairing, e.g. GC→GU). |  |  |  | dna, rna, codon |

## RNA base pairs (per pair)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `pair_i` | 5' column of the base pair (1-based). |  |  |  | all |
| `pair_j` | 3' column of the base pair (1-based). |  |  |  | all |
| `sequences_complete` | Sequences with a nucleotide at both columns. |  |  |  | all |
| `watson_crick` | Sequences with AU, UA, GC or CG. |  |  |  | all |
| `wobble` | Sequences with GU or UG. |  |  |  | all |
| `noncanonical` | Sequences with any other nucleotide pair. |  |  |  | all |
| `incomplete` | Sequences with a gap or ambiguity at either column. |  |  |  | all |
| `canonical_fraction` | (watson_crick + wobble) / sequences_complete. |  |  | No complete sequence. | all |
| `dominant_pair` | Most common canonical pair; reference for covariation. |  |  |  | all |
| `compensatory_sequences` | Canonical pairs differing from dominant_pair at both positions. |  |  |  | all |
| `consistent_sequences` | Canonical pairs differing from dominant_pair at one position. |  |  |  | all |
| `covarying` | True when compensatory_sequences > 0. |  |  |  | all |

## Profile match/insert states (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `profile_source` | Where match/insert states came from (A2M/A3M or #=GC RF). |  |  |  | all |
| `match_states` | Match-state columns. |  |  |  | all |
| `insert_columns` | Insert columns (#=GC RF '.'). |  |  |  | all |
| `mean_match_state_coverage` | Mean fraction of match states occupied by a residue. |  |  |  | all |
| `total_inserted_residues` | Residues in insert states, all sequences. |  |  |  | all |
| `sequences_with_insertions` | Sequences with at least one inserted residue. |  |  |  | all |
| `insert_positions_used` | Positions between match states where any sequence inserts. |  |  |  | all |
| `max_sequences_inserting_at_one_position` | Largest number of sequences inserting at one position. |  |  |  | all |

## Profile match/insert states (per sequence)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_id` | Sequence identifier. |  |  |  | all |
| `match_states` | Match-state columns. |  |  |  | all |
| `match_states_occupied` | Match states with a residue in this sequence. |  |  |  | all |
| `match_state_coverage` | match_states_occupied / match_states. |  |  |  | all |
| `deletions` | Match states deleted ('-') in this sequence. |  |  |  | all |
| `inserted_residues` | Lowercase (insert-state) residues in this sequence. |  |  |  | all |
| `insertion_events` | Positions with at least one inserted residue. |  |  |  | all |

## Partitions (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `partitions` | Declared partitions. |  |  |  | all |
| `partition_sites_assigned` | Alignment columns in at least one partition. |  |  |  | all |
| `partition_sites_unassigned` | Columns in no partition. |  |  |  | all |
| `partition_sites_overlapping` | Column assignments beyond the first (overlapping partitions, e.g. genes and codon positions). |  |  |  | all |
| `partitions_failing_composition_test` | Partitions with composition χ² p < 0.05. |  |  |  | all |

## Partitions (per partition)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `partition` | Partition name. |  |  |  | all |
| `datatype` | Data type from the partition file, if given. |  |  |  | all |
| `sites_declared` | Columns listed for the partition. |  |  |  | all |
| `sites_in_alignment` | Declared columns inside the alignment. |  |  |  | all |
| `occupancy` | Residue fraction of the partition's cells. |  |  |  | all |
| `gap_fraction` | Gap fraction of the partition's cells. |  |  |  | all |
| `sequences_without_data` | Sequences with only gaps/missing data in the partition. |  |  |  | all |
| `variable_sites` | Variable columns in the partition. |  |  |  | all |
| `parsimony_informative_sites` | Parsimony-informative columns in the partition. |  |  |  | all |
| `entropy_mean` | Mean column entropy in the partition. |  |  |  | all |
| `gc_content` | GC fraction (nucleotides). |  |  |  | dna, rna, codon |
| `composition_chi2` | Pearson χ² of the sequence × residue count table. |  |  | < 2 sequences or states with data. | all |
| `composition_chi2_df` | Degrees of freedom (rows − 1)(columns − 1). |  |  |  | all |
| `composition_chi2_p` | Upper-tail χ² p-value; small values indicate compositional heterogeneity. |  |  |  | all |

## Block alignments (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `alignment_kind` | Data model of the input. |  |  |  | all |
| `detected_format` | Detected input format. |  |  |  | all |
| `blocks` | Alignment blocks. |  |  |  | all |
| `genomes` | Distinct genomes (MAF 'genome.chrom' prefix; XMFA sequence index/file). |  |  |  | all |
| `reference_genome` | Genome of the first row of the first block. |  |  |  | all |
| `total_block_columns` | Sum of block widths. |  |  |  | all |
| `min_block_columns` | Narrowest block. |  |  |  | all |
| `median_block_columns` | Median block width. |  |  |  | all |
| `max_block_columns` | Widest block. |  |  |  | all |
| `n50_block_columns` | Block-width N50. |  |  |  | all |
| `mean_genomes_per_block` | Mean distinct genomes per block. |  |  |  | all |
| `blocks_with_all_genomes` | Blocks containing every genome. |  |  |  | all |
| `reference_aligned_fraction` | Reference bases covered by blocks / reference size. |  |  | Source sizes not in the format (XMFA, AXT). | all |
| `minus_strand_rows` | Block rows on the minus strand. |  |  |  | all |
| `collinearity_breakpoints` | Adjacent blocks (reference order) whose chromosome, strand or order differs in another genome. |  |  |  | all |
| `duplicated_bases` | Aligned bases minus covered bases (bases aligned more than once). |  |  |  | all |

## Block alignments (per block)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `block` | Block number (1-based, file order). |  |  |  | all |
| `columns` | Block width. |  |  |  | all |
| `rows` | Rows in the block. |  |  |  | all |
| `genomes` | Distinct genomes in the block. |  |  |  | all |
| `reference_src` | Source sequence of the first row. |  |  |  | all |
| `reference_start` | 0-based forward-strand start of the first row. |  |  |  | all |
| `reference_end` | Half-open end of the first row. |  |  |  | all |
| `reference_strand` | Strand of the first row. |  |  |  | all |
| `gap_fraction` | Gap fraction of the block. |  |  |  | all |
| `variable_sites` | Variable columns in the block. |  |  |  | all |
| `mean_consensus_identity` | Mean identity of rows to the block consensus. |  |  |  | all |
| `score` | Block score from the file. |  |  |  | all |

## Block alignments (per genome)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `genome` | Genome name. |  |  |  | all |
| `is_reference` | True for the reference genome. |  |  |  | all |
| `blocks_present` | Blocks containing this genome. |  |  |  | all |
| `block_fraction` | blocks_present / blocks. |  |  |  | all |
| `aligned_bases` | Sum of aligned interval lengths. |  |  |  | all |
| `covered_bases` | Union of aligned intervals. |  |  |  | all |
| `duplicated_bases` | aligned_bases − covered_bases. |  |  |  | all |
| `genome_size` | Sum of source sizes (MAF). |  |  |  | all |
| `aligned_fraction` | covered_bases / genome_size. |  |  |  | all |
| `minus_strand_rows` | Rows on the minus strand. |  |  |  | all |
| `blocks_with_multiple_copies` | Blocks containing this genome more than once (paralogy/duplication). |  |  |  | all |
| `collinearity_breakpoints` | Rearrangement breakpoints relative to the reference. |  |  |  | all |

## Read alignments (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `alignment_kind` | Data model of the input. |  |  |  | all |
| `detected_format` | Detected input format. |  |  |  | all |
| `references` | Reference sequences in the header. |  |  |  | all |
| `reference_bases` | Total reference length. |  |  |  | all |
| `records` | Alignment records. |  |  |  | all |
| `primary_reads` | Records that are neither secondary nor supplementary. |  |  |  | all |
| `secondary_records` | Flag 0x100. |  |  |  | all |
| `supplementary_records` | Flag 0x800. |  |  |  | all |
| `unmapped_reads` | Primary records with flag 0x4. |  |  |  | all |
| `mapping_rate` | Mapped primary reads / primary reads. |  |  |  | all |
| `duplicate_reads` | Primary reads with flag 0x400. |  |  |  | all |
| `qc_fail_reads` | Primary reads with flag 0x200. |  |  |  | all |
| `paired_reads` | Primary reads with flag 0x1. |  |  |  | all |
| `properly_paired_reads` | Primary reads with flag 0x2. |  |  |  | all |
| `mean_mapq` | Mean MAPQ of mapped primary reads (255 excluded). |  |  |  | all |
| `median_mapq` | Median MAPQ. |  |  |  | all |
| `mapq_ge_30_fraction` | Fraction with MAPQ ≥ 30. |  |  |  | all |
| `mapq_0_fraction` | Fraction with MAPQ 0 (multi-mapping). |  |  |  | all |
| `aligned_bases` | CIGAR M/=/X bases of mapped primary reads. |  |  |  | all |
| `inserted_bases` | CIGAR I bases. |  |  |  | all |
| `deleted_bases` | CIGAR D bases. |  |  |  | all |
| `indel_rate` | (I + D) / aligned_bases. |  |  |  | all |
| `mismatch_rate` | Mismatches / aligned bases, from =/X CIGAR or NM − I − D. |  |  | NM missing and no =/X CIGAR. | all |
| `soft_clipped_bases` | CIGAR S bases. |  |  |  | all |
| `hard_clipped_bases` | CIGAR H bases. |  |  |  | all |
| `clipped_read_fraction` | Mapped primary reads with any clipping. |  |  |  | all |
| `spliced_reads` | Mapped reads with CIGAR N (spliced/split). |  |  |  | all |
| `mean_read_length` | Mean read length. |  |  |  | all |
| `n50_read_length` | Read-length N50. |  |  |  | all |
| `mean_base_quality` | Mean Phred base quality. |  |  |  | all |
| `median_insert_size` | Median \|TLEN\| of properly paired first mates. |  |  |  | all |
| `mean_depth` | Length-weighted mean depth over references. |  |  | No @SQ lengths, or references too large. | all |
| `breadth_1x` | Fraction of reference bases with depth ≥ 1. |  |  |  | all |
| `breadth_10x` | Fraction of reference bases with depth ≥ 10. |  |  |  | all |

## Read alignments (per reference)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `reference` | Reference name. |  |  |  | all |
| `length` | Reference length (@SQ LN). |  |  |  | all |
| `primary_mapped_reads` | Primary reads mapped to the reference. |  |  |  | all |
| `mean_depth` | Mean depth. |  |  |  | all |
| `median_depth` | Median depth. |  |  |  | all |
| `breadth_1x` | Fraction of bases with depth ≥ 1. |  |  |  | all |
| `breadth_10x` | Fraction of bases with depth ≥ 10. |  |  |  | all |
| `breadth_30x` | Fraction of bases with depth ≥ 30. |  |  |  | all |

## Mappings and hits (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `alignment_kind` | Data model of the input. |  |  |  | all |
| `detected_format` | Detected input format. |  |  |  | all |
| `records` | Hit/mapping records. |  |  |  | all |
| `queries` | Distinct queries. |  |  |  | all |
| `targets` | Distinct targets (sequences or graph paths). |  |  |  | all |
| `multi_hit_queries` | Queries with more than one record. |  |  |  | all |
| `mean_identity` | Mean per-record identity (format-specific definition). |  |  |  | all |
| `median_identity` | Median per-record identity. |  |  |  | all |
| `length_weighted_identity` | Identity weighted by alignment length. |  |  |  | all |
| `total_aligned_length` | Sum of alignment lengths. |  |  |  | all |
| `median_alignment_length` | Median alignment length. |  |  |  | all |
| `minus_strand_fraction` | Records on the minus strand. |  |  |  | all |
| `mean_mapq` | Mean MAPQ (PAF/GAF). |  |  |  | all |
| `mapq_ge_30_fraction` | Fraction with MAPQ ≥ 30. |  |  |  | all |
| `mean_query_coverage` | Mean over queries of covered bases / query length. |  |  | Query lengths not in the format. | all |
| `mean_target_coverage` | Mean over targets of covered bases / target length. |  |  | Target lengths not in the format. | all |
| `min_evalue` | Smallest e-value (BLAST). |  |  |  | all |

## Mappings and hits (per query / target)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `query` | Query name. |  |  |  | all |
| `target` | Target name. |  |  |  | all |
| `length` | Declared length. |  |  |  | all |
| `hits` | Records involving this sequence. |  |  |  | all |
| `covered_bases` | Union of aligned intervals. |  |  |  | all |
| `coverage` | covered_bases / length. |  |  |  | all |

## Sequence graphs (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `alignment_kind` | Data model of the input. |  |  |  | all |
| `detected_format` | Detected input format. |  |  |  | all |
| `segments` | Segments (S lines). |  |  |  | all |
| `total_segment_bases` | Sum of segment lengths. |  |  |  | all |
| `mean_segment_length` | Mean segment length. |  |  |  | all |
| `n50_segment_length` | Segment-length N50. |  |  |  | all |
| `links` | Links (L lines). |  |  |  | all |
| `connected_components` | Weakly connected components. |  |  |  | all |
| `isolated_segments` | Segments without links. |  |  |  | all |
| `tip_segments` | Segments with links on one side only. |  |  |  | all |
| `branching_segment_ends` | Segment ends with more than one link (bubble/branch points). |  |  |  | all |
| `max_degree` | Largest number of links on one segment. |  |  |  | all |
| `paths` | Paths (P) and walks (W). |  |  |  | all |
| `mean_path_bases` | Mean path length in bases. |  |  |  | all |
| `segments_on_paths` | Segments traversed by at least one path. |  |  |  | all |
| `core_segments` | Segments on every path. |  |  |  | all |
| `shell_segments` | Segments on more than one but not all paths. |  |  |  | all |
| `private_segments` | Segments on exactly one path. |  |  |  | all |
| `core_bases_fraction` | Bases in core segments / all segment bases. |  |  |  | all |
| `reference_rank0_segments` | rGFA segments with SR:i:0 (reference backbone). |  |  |  | all |

## Sequence graphs (per segment)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `segment` | Segment name. |  |  |  | all |
| `length` | Segment length. |  |  |  | all |
| `degree` | Links touching the segment. |  |  |  | all |
| `left_links` | Links on the segment's left end. |  |  |  | all |
| `right_links` | Links on the segment's right end. |  |  |  | all |
| `paths` | Paths/walks traversing the segment. |  |  |  | all |
| `rgfa_rank` | rGFA SR rank. |  |  |  | all |

## Sequence graphs (per path)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `path` | Path or walk name (walk: sample#hap#seq:start-end). |  |  |  | all |
| `segments` | Segment steps. |  |  |  | all |
| `bases` | Path length in bases. |  |  |  | all |
| `distinct_segments` | Distinct segments visited. |  |  |  | all |

## Profile HMMs (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `alignment_kind` | Data model of the input. |  |  |  | all |
| `detected_format` | Detected input format. |  |  |  | all |
| `model_name` | NAME line. |  |  |  | all |
| `model_accession` | ACC line. |  |  |  | all |
| `alphabet` | ALPH line. |  |  |  | all |
| `model_length` | Match states (LENG). |  |  |  | all |
| `training_sequences` | Sequences used to build the model (NSEQ). |  |  |  | all |
| `effective_sequences` | Effective sequence number (EFFN / NEFF). |  |  |  | all |
| `mean_match_entropy_bits` | Mean Shannon entropy of match-state emissions (HMMER3). |  |  |  | all |
| `mean_match_relative_entropy_bits` | Mean relative entropy of match emissions vs the model COMPO composition. |  |  |  | all |

## Profile HMMs (per match state)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `match_state` | Match-state index. |  |  |  | all |
| `consensus` | Most probable residue. |  |  |  | all |
| `max_emission` | Largest emission probability. |  |  |  | all |
| `entropy_bits` | Shannon entropy of emissions. |  |  |  | all |
| `relative_entropy_bits` | Relative entropy vs COMPO. |  |  |  | all |

## Completeness (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `alignment_completeness` | Fraction of all cells observed (AliStat Ca). *Observed = canonical residue; gaps, missing, ambiguity codes and padding are unobserved.* | `observed cells / (N L)` | 0-1 |  | all |
| `min_sequence_completeness` | Lowest per-sequence completeness (AliStat Cr). |  | 0-1 |  | all |
| `min_site_completeness` | Lowest per-column completeness (AliStat Cc). |  | 0-1 |  | all |
| `sequence_completeness_quantiles` | Per-sequence completeness. Quantiles q05, q25, q50, q75, q95. |  |  |  | all |
| `site_completeness_quantiles` | Per-column completeness. Quantiles q05, q25, q50, q75, q95. |  |  |  | all |
| `min_pair_completeness_threshold` | Threshold used for pairs_below_threshold. |  |  |  | all |
| `pair_completeness_method` | exact, or why pairwise overlap was not computed. |  |  |  | all |
| `mean_pair_completeness` | Mean AliStat C_ij over pairs. Cost: O(N^2 L) time, O(N^2) memory; skipped above max_exact_sequences. | `C_ij = columns observed in both i and j / L` |  |  | all |
| `min_pair_completeness` | Smallest C_ij: the least-overlapping pair. |  |  |  | all |
| `pair_completeness_quantiles` | C_ij over all pairs. Quantiles q05, q25, q50, q75, q95. |  |  |  | all |
| `pairs_below_threshold` | Pairs with C_ij below the threshold. |  |  |  | all |
| `pairs_without_overlap` | Pairs sharing no observed column (distance undefined). |  |  |  | all |
| `fraction_pairs_below_threshold` | pairs_below_threshold / pairs. |  |  |  | all |
| `min_shared_sites` | Fewest columns observed in both sequences of any pair. |  |  |  | all |

## Completeness (per sequence)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_id` | Sequence identifier. |  |  |  | all |
| `sequence_completeness` | Observed columns / L (AliStat Cr). |  |  |  | all |
| `least_overlapping_partner` | Sequence with the smallest C_ij to this one. |  |  |  | all |
| `min_pair_completeness` | That smallest C_ij. |  |  |  | all |
| `pairs_below_threshold` | Partners with C_ij below the threshold. |  |  |  | all |

## Gap structure (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `gap_openings` | Internal gap runs (a run starts after a non-gap cell). *Terminal gaps (before the first / after the last non-gap cell) are counted separately, not as openings.* |  |  |  | all |
| `gap_extensions` | Internal gap cells beyond the first of each run. | `internal_gap_cells - gap_openings` |  |  | all |
| `internal_gap_cells` | Gap cells between a sequence's first and last non-gap cell. |  |  |  | all |
| `terminal_gap_cells` | Gap cells before the first / after the last non-gap cell. |  |  |  | all |
| `terminal_gap_fraction` | terminal_gap_cells / all gap cells. |  |  | no gaps | all |
| `mean_gap_run_length` | Mean internal gap-run length. |  |  | no internal gaps | all |
| `median_gap_run_length` | Median internal gap-run length. |  |  |  | all |
| `max_gap_run_length` | Longest internal gap run. |  |  |  | all |
| `single_site_gap_runs` | Internal gap runs of length 1. |  |  |  | all |
| `longest_gap_free_block` | Longest run of columns containing no gap. |  |  |  | all |
| `gap_free_blocks` | Number of maximal gap-free column runs. |  |  |  | all |
| `longest_complete_block` | Longest run of columns observed in every sequence. |  |  |  | all |
| `complete_blocks` | Number of maximal complete column runs. |  |  |  | all |
| `distinct_column_gap_patterns` | Distinct gap-presence patterns among gapped columns. |  |  |  | all |
| `distinct_sequence_gap_patterns` | Distinct gap-position sets among sequences. |  |  |  | all |
| `mean_pairwise_gap_concordance` | Mean Jaccard similarity of gap positions over pairs. | `\|G_i ∩ G_j\| / \|G_i ∪ G_j\| (1 when neither has gaps)` | 0-1 |  | all |
| `fragmentation_index` | Mean over sequences of 1 - longest observed segment / observed cells; 0 when every sequence is one contiguous block. |  | 0-1 |  | all |

## Gap structure (per sequence)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_id` | Sequence identifier. |  |  |  | all |
| `internal_gap_openings` | Internal gap runs. |  |  |  | all |
| `internal_gap_extensions` | Internal gap cells beyond each run's first. |  |  |  | all |
| `internal_gap_cells` | Internal gap cells. |  |  |  | all |
| `terminal_gap_cells` | Terminal gap cells. |  |  |  | all |
| `longest_internal_gap` | Longest internal gap run. |  |  |  | all |
| `observed_segments` | Maximal runs of observed cells. |  |  |  | all |
| `longest_observed_segment` | Longest such run. |  |  |  | all |
| `fragmentation` | 1 - longest_observed_segment / observed cells. |  |  |  | all |
| `gap_pattern_distance` | 1 - mean Jaccard similarity of this sequence's gap positions to every other sequence's. |  |  |  | all |
| `gap_pattern_robust_z` | Robust z of gap_pattern_distance (median/MAD); large positive values mark unusual gap placement. |  |  | MAD = 0 | all |

## Gap-run length distribution

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `run_length` | Internal gap-run length. |  |  |  | all |
| `runs` | Runs of that length. |  |  |  | all |

## Contiguous blocks

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `block_type` | gap_free (no gap in the column) or complete (all observed). |  |  |  | all |
| `start` | First column (1-based). |  |  |  | all |
| `end` | Last column (1-based, inclusive). |  |  |  | all |
| `length` | Columns in the block. |  |  |  | all |

## Redundancy and N_eff (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `identity_threshold` | Identity at or above which two sequences are neighbours. |  |  |  | all |
| `min_overlap` | Neighbours' comparable sites must cover at least this fraction of the longer sequence's residues. |  |  |  | all |
| `effective_sequences` | Effective number of sequences N_eff. Assumes: Identity over sites where both sequences have a residue. Cost: O(N^2 L) time, O(N^2) memory; skipped above max_exact_sequences. | `sum_i 1/n_i, n_i = neighbours of i (including i)` | 1-N |  | all |
| `effective_sequence_fraction` | N_eff / N. |  | 0-1 |  | all |
| `effective_sequences_per_site` | N_eff / L. |  |  |  | all |
| `mean_neighbours` | Mean n_i. |  |  |  | all |
| `max_neighbours` | Largest n_i. |  |  |  | all |
| `clusters` | Connected components of the neighbour graph (single linkage). |  |  |  | all |
| `largest_cluster` | Size of the largest cluster. |  |  |  | all |
| `singleton_clusters` | Clusters with one sequence. |  |  |  | all |
| `pairwise_identity_quantiles` | Pairwise identity. Quantiles q05, q25, q50, q75, q95. |  |  |  | all |

## Sequence weights and clusters

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_id` | Sequence identifier. |  |  |  | all |
| `neighbours` | n_i. |  |  |  | all |
| `weight` | w_i = 1/n_i; weights sum to N_eff. |  |  |  | all |
| `cluster` | Cluster label (0-based, first-seen order). |  |  |  | all |
| `cluster_size` | Sequences in the cluster. |  |  |  | all |

## Cluster size distribution

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `cluster_size` | Cluster size. |  |  |  | all |
| `clusters` | Clusters of that size. |  |  |  | all |

## Site diversity and patterns (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `distinct_site_patterns` | Distinct columns (all characters, gaps included). |  |  |  | all |
| `site_pattern_ratio` | distinct_site_patterns / L. |  | 0-1 |  | all |
| `repeated_pattern_sites` | L - distinct_site_patterns. |  |  |  | all |
| `max_pattern_multiplicity` | Most frequent column pattern's count. |  |  |  | all |
| `site_pattern_entropy` | Shannon entropy (bits) of the pattern frequency distribution. |  |  |  | all |
| `mean_gini_simpson` | Mean over observed columns of 1 - sum p_c^2. *Canonical residues only.* |  |  |  | all |
| `mean_effective_alphabet_size` | Mean 2^H over observed columns. |  |  |  | all |
| `mean_weighted_entropy` | Mean column entropy with sequence weights w_i. |  |  | no exact pairwise identities | all |
| `weighted` | Whether sequence weights were available. |  |  |  | all |
| `sites_1_state` | Columns with one observed state. |  |  |  | all |
| `sites_2_states` | Columns with two observed states. |  |  |  | all |
| `sites_3_states` | Columns with three observed states. |  |  |  | all |
| `sites_4plus_states` | Columns with four or more observed states. |  |  |  | all |
| `singleton_minor_sites` | Variable columns whose non-majority states total 1. |  |  |  | all |
| `doubleton_minor_sites` | Variable columns whose non-majority states total 2. |  |  |  | all |
| `minor_count_3plus_sites` | Variable columns whose non-majority states total >= 3. |  |  |  | all |

## Site diversity (per site)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `site_index` | Column (1-based). |  |  |  | all |
| `simpson` | sum p_c^2. |  |  |  | all |
| `gini_simpson` | 1 - sum p_c^2. |  |  |  | all |
| `effective_alphabet_size` | 2^H. |  |  |  | all |
| `weighted_entropy` | Entropy (bits) with sequence weights. |  |  |  | all |
| `minor_state_count` | Observed residues not in the majority state. |  |  |  | all |
| `pattern_multiplicity` | Columns sharing this column's pattern. |  |  |  | all |
| `pattern_first_site` | First column with this pattern. |  |  |  | all |

## Composition heterogeneity (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `rcv` | Relative composition variability (Phillips & Penny 2003). *Canonical residues only.* | `sum_c sum_i \|x_ic - mean_c\| / (N L)` |  |  | all |
| `max_rcvt` | Largest per-sequence RCV contribution. |  |  |  | all |
| `sequences_tested` | Sequences with at least one residue. |  |  |  | all |
| `sequences_failing_composition_test` | Sequences with χ² p < alpha. |  |  |  | all |
| `fraction_failing_composition_test` | Failing / tested. |  |  |  | all |
| `alpha` | Significance level. |  |  |  | all |
| `mean_composition_jsd` | Mean Jensen–Shannon divergence (bits) of sequence composition from the alignment composition. |  | 0-1 |  | all |
| `max_composition_jsd` | Largest such divergence. |  |  |  | all |

## Composition tests (per sequence)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_id` | Sequence identifier. |  |  |  | all |
| `rcvt` | Per-sequence RCV contribution: sum_c \|x_ic - mean_c\| / (N L). |  |  |  | all |
| `composition_chi2` | Pearson χ² of residue counts vs alignment frequencies. Assumes: Sites treated as independent; a screening test, as in IQ-TREE. |  |  |  | all |
| `composition_df` | Degrees of freedom: observed states - 1. |  |  |  | all |
| `composition_p` | χ² p-value. |  |  |  | all |
| `composition_test_failed` | p < alpha. |  |  |  | all |
| `composition_jsd` | Jensen–Shannon divergence (bits) from the alignment composition. |  |  |  | all |

## Symmetry tests (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `pair_method` | exact (all pairs) or sampled. |  |  |  | all |
| `pairs_tested` | Pairs tested. |  |  |  | all |
| `pairs_total` | N(N-1)/2. |  |  |  | all |
| `sample_seed` | Random seed when sampled. |  |  |  | all |
| `alpha` | Significance level. |  |  |  | all |
| `fraction_pairs_rejecting_bowker` | Pairs with Bowker p < alpha (rejecting symmetry: stationary, reversible, homogeneous evolution). Assumes: No multiple-testing correction; see Naser-Khdour et al. 2019. |  |  |  | all |
| `fraction_pairs_rejecting_marginal` | Pairs with Stuart–Maxwell marginal symmetry p < alpha (rejecting stationarity). |  |  |  | all |
| `fraction_pairs_rejecting_internal` | Pairs with internal symmetry p < alpha. |  |  |  | all |
| `min_bowker_p` | Smallest Bowker p-value. |  |  |  | all |
| `multiple_testing` | Multiple-testing note. |  |  |  | all |

## Symmetry tests (per pair)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_a` | First sequence. |  |  |  | all |
| `sequence_b` | Second sequence. |  |  |  | all |
| `sites` | Sites where both carry a canonical residue (divergence-matrix total). |  |  |  | all |
| `bowker` | Bowker statistic sum_{i<j} (n_ij - n_ji)^2 / (n_ij + n_ji). |  |  |  | all |
| `bowker_df` | Off-diagonal pairs with n_ij + n_ji > 0. |  |  |  | all |
| `bowker_p` | χ² p-value. |  |  |  | all |
| `marginal` | Stuart–Maxwell statistic d' V^-1 d on the first k-1 margins. |  |  |  | all |
| `marginal_df` | Rank of V. |  |  |  | all |
| `marginal_p` | χ² p-value. |  |  |  | all |
| `internal` | bowker - marginal. |  |  |  | all |
| `internal_df` | bowker_df - marginal_df. |  |  |  | all |
| `internal_p` | χ² p-value. |  |  | internal_df <= 0 | all |

## Windowed statistics (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `window` | Window width (sites; a multiple of 3 for codon alignments). |  |  |  | all |
| `step` | Step between window starts. |  |  |  | all |
| `windows` | Number of windows. |  |  |  | all |
| `most_variable_window` | Window with the highest variable_fraction. |  |  |  | all |
| `gappiest_window` | Window with the highest gap_fraction. |  |  |  | all |
| `highest_entropy_window` | Window with the highest mean_entropy. |  |  |  | all |
| `windows_with_outlier_sequences` | Windows with at least one outlier sequence. |  |  |  | all |
| `outlier_z_threshold` | Robust-z threshold for window outliers. |  |  |  | all |
| `outlier_mad_floor` | Lower bound on the scaled MAD (identity units). |  |  |  | all |

## Protein properties (per site)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `site_index` | Column (1-based). |  |  |  | all |
| `hydropathy_sd` | SD of Kyte–Doolittle hydropathy over residues in the column. |  |  |  | all |
| `volume_sd` | SD of residue volume (Zamyatnin 1972, Å^3). |  |  |  | all |
| `charge_variance` | Variance of net charge (K, R +1; D, E -1; H +0.1). |  |  |  | all |
| `polar_fraction` | Fraction of residues that are polar (STNQCYHKRDE). |  |  |  | all |
| `substitutions_vs_consensus` | Residues differing from the column consensus. |  |  |  | all |
| `conservative_substitutions` | Of those, BLOSUM62(residue, consensus) > 0. |  |  |  | all |
| `radical_substitutions` | Of those, BLOSUM62(residue, consensus) < 0. |  |  |  | all |
| `class_changes` | Of those, a change of residue class (hydrophobic, polar, positive, negative, special G/P/C). |  |  |  | all |
| `radical_fraction` | radical / substitutions. |  |  |  | all |

## Protein substitution scores (per sequence)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_id` | Sequence identifier. |  |  |  | all |
| `mean_blosum62_to_consensus` | Mean BLOSUM62 score of residues vs the column consensus. |  |  |  | all |
| `substitution_robust_z` | Robust z of that mean; strongly negative = unusually radical replacements. |  |  |  | all |
| `radical_substitutions` | Radical substitutions in the sequence. |  |  |  | all |
| `conservative_substitutions` | Conservative substitutions in the sequence. |  |  |  | all |

## Protein properties (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `substitutions_vs_consensus` | Residues differing from their column consensus. |  |  |  | all |
| `conservative_fraction` | BLOSUM62 > 0 share of those substitutions. |  |  |  | all |
| `radical_fraction` | BLOSUM62 < 0 share. |  |  |  | all |
| `neutral_fraction` | BLOSUM62 = 0 share. |  |  |  | all |
| `class_change_fraction` | Share changing residue class. |  |  |  | all |
| `median_hydropathy_sd` | Median column hydropathy SD. |  |  |  | all |
| `sequences_with_substitution_z_below_minus_3_5` | Sequences with substitution robust z < -3.5. |  |  |  | all |

## Codon evolution (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `genetic_code` | NCBI genetic code table. |  |  |  | all |
| `codon_sites` | Codon columns (L / 3). |  |  |  | all |
| `variable_codon_sites` | Codon columns with more than one complete, non-stop codon. |  |  |  | all |
| `synonymous_only_codon_sites` | Variable codon columns whose codons all encode one amino acid. |  |  |  | all |
| `nonsynonymous_codon_sites` | Variable codon columns with more than one amino acid. |  |  |  | all |
| `effective_number_of_codons` | Wright's ENC from all complete codons. |  | 20-61 | a synonymous family class is unobserved | all |
| `fourfold_sites` | Third positions whose consensus codon is four-fold degenerate. |  |  |  | all |
| `fourfold_nucleotide_diversity` | Mean pairwise difference per four-fold site among codons sharing the consensus first two bases. |  |  |  | all |
| `dnds_method` | Estimator used for dN and dS. |  |  |  | all |
| `dnds_pair_method` | exact (all pairs) or sampled. |  |  |  | all |
| `dnds_pairs` | Pairs evaluated. |  |  |  | all |
| `dnds_pairs_defined` | Pairs with a defined dN/dS. |  |  |  | all |
| `median_pairwise_dn_ds` | Median pairwise dN/dS. Assumes: Counting estimate between pairs; not a phylogenetic selection test. |  |  |  | all |
| `mean_pairwise_ds` | Mean pairwise dS. |  |  |  | all |
| `mean_pairwise_dn` | Mean pairwise dN. |  |  |  | all |
| `dnds_pairs_saturated` | Pairs with p >= 0.75 (Jukes–Cantor undefined). |  |  |  | all |

## Pairwise dN/dS

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_a` | First sequence. |  |  |  | all |
| `sequence_b` | Second sequence. |  |  |  | all |
| `codons` | Codons compared (complete, non-stop in both). |  |  |  | all |
| `synonymous_sites` | S (Nei–Gojobori sites, averaged over the two sequences). |  |  |  | all |
| `nonsynonymous_sites` | N = 3 codons - S. |  |  |  | all |
| `synonymous_differences` | S_d. |  |  |  | all |
| `nonsynonymous_differences` | N_d. |  |  |  | all |
| `ps` | S_d / S. |  |  |  | all |
| `pn` | N_d / N. |  |  |  | all |
| `ds` | Jukes–Cantor corrected p_S. | `-3/4 ln(1 - 4 p_S / 3)` |  |  | all |
| `dn` | Jukes–Cantor corrected p_N. |  |  |  | all |
| `dn_ds` | dN / dS. |  |  | dS = 0, saturation, or < 10 codons | all |
| `note` | Why a value is undefined. |  |  |  | all |

## Coevolution (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `candidate_sites` | Columns analysed: variable, occupancy >= min_occupancy, the max_sites most entropic. |  |  |  | all |
| `max_sites` | Cap on candidate columns. Cost: O(N (S k)^2) time. |  |  |  | all |
| `min_occupancy` | Minimum fraction of sequences observed in a candidate column. |  |  |  | all |
| `pseudocount` | lambda in p = (1 - lambda) f + lambda / q. |  |  |  | all |
| `weighting` | Sequence weighting used. |  |  |  | all |
| `effective_sequences` | Sum of weights. |  |  |  | all |
| `min_separation` | Minimum column distance for reported pairs. |  |  |  | all |
| `pairs` | Column pairs evaluated. |  |  |  | all |
| `note` | Why no pairs were evaluated. |  |  |  | all |
| `mean_mi` | Mean MI over evaluated pairs (bits). |  |  |  | all |
| `max_mi_apc` | Largest APC-corrected MI. |  |  |  | all |
| `top_pair` | Columns of the highest-ranked pair. |  |  |  | all |
| `pairs_reported` | Rows in coevolution_pairs.csv. |  |  |  | all |

## Coevolution (top column pairs)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `rank` | Rank by APC-corrected MI. |  |  |  | all |
| `site_i` | First column (1-based). |  |  |  | all |
| `site_j` | Second column (1-based). |  |  |  | all |
| `mi` | Mutual information (bits). | `sum p(x,y) log2 p(x,y) / (p(x) p(y))` |  |  | all |
| `mi_apc` | MI - MI(i,.) MI(.,j) / MI(.,.) (Dunn et al. 2008). |  |  |  | all |
| `normalized_mi` | 2 MI / (H_i + H_j). |  | 0-1 |  | all |
| `joint_entropy` | H(i, j) in bits. |  |  |  | all |
| `conditional_entropy_j_given_i` | H(j \| i) = H(i, j) - H(i). |  |  |  | all |
| `effective_pair_weight` | Weighted sequences with residues in both columns. |  |  |  | all |

## Alignment comparison (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `pair_method` | exact (all sequence pairs) or sampled. |  |  |  | all |
| `sequence_pairs` | Sequence pairs compared. |  |  |  | all |
| `reference_residue_pairs` | Aligned residue pairs implied by the reference. |  |  |  | all |
| `test_residue_pairs` | Aligned residue pairs implied by the test alignment. |  |  |  | all |
| `shared_residue_pairs` | Pairs implied by both. |  |  |  | all |
| `sp_score` | Sum-of-pairs score: shared / reference pairs (recall). |  | 0-1 |  | all |
| `precision` | shared / test pairs. |  | 0-1 |  | all |
| `f1` | Harmonic mean of SP and precision. |  |  |  | all |
| `tc_score` | Total-column score: reference columns (>= 2 residues) reproduced exactly. |  | 0-1 |  | all |
| `reference_columns` | Reference columns with at least two residues. |  |  |  | all |
| `mean_column_agreement` | Mean column agreement over test columns. |  |  |  | all |
| `columns_below_0_5_agreement` | Test columns with agreement < 0.5. |  |  |  | all |

## Alignment comparison (per column)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `site_index` | Test-alignment column (1-based). |  |  |  | all |
| `residue_pairs` | Residue pairs in the column (over compared sequence pairs). |  |  |  | all |
| `agreement` | Fraction of those pairs also aligned in the other alignment. |  |  |  | all |

## Alignment comparison (per sequence)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `sequence_id` | Sequence identifier. |  |  |  | all |
| `pairs_recovered` | Fraction of the reference's residue pairs involving this sequence that the test alignment recovers. |  |  |  | all |

## Alignment ensemble (summary)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `alternatives` | Alternative alignments compared. |  |  |  | all |
| `mean_sp_vs_alternatives` | Mean SP of the test alignment against each alternative. |  |  |  | all |
| `mean_tc_vs_alternatives` | Mean TC against each alternative. |  |  |  | all |
| `mean_column_confidence` | Mean per-column confidence. |  |  |  | all |
| `columns_below_0_5_confidence` | Columns with confidence < 0.5. |  |  |  | all |
| `columns_fully_consistent` | Columns aligned identically in every alternative. |  |  |  | all |

## Alignment ensemble (per column)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `site_index` | Test-alignment column (1-based). |  |  |  | all |
| `confidence` | Mean column agreement over the alternatives. |  |  |  | all |
| `min_agreement` | Lowest agreement with any alternative. |  |  |  | all |

## Windowed statistics (per window)

| Field | Definition | Formula | Range | Undefined when | Applies to |
|---|---|---|---|---|---|
| `window_start` | First column (1-based). |  |  |  | all |
| `window_end` | Last column (inclusive). |  |  |  | all |
| `sites` | Columns in the window. |  |  |  | all |
| `occupancy` | Mean column occupancy. |  |  |  | all |
| `gap_fraction` | Mean column gap fraction. |  |  |  | all |
| `variable_fraction` | Fraction of variable columns. |  |  |  | all |
| `parsimony_informative_fraction` | Fraction of parsimony-informative columns. |  |  |  | all |
| `mean_entropy` | Mean column entropy (bits). |  |  |  | all |
| `mean_consensus_identity` | Mean identity to the consensus in the window, over sequences with at least half the window observed. |  |  |  | all |
| `outlier_sequences` | Sequences with window consensus-identity robust z < threshold. |  |  |  | all |
| `gc_content` | GC fraction (nucleotides). |  |  |  | dna, rna, codon |

