"""Metric definitions for the characterization analyses and alignment comparison.

Levels are the analysis names (summaries) and table names (rows). Merged into the
registry by :mod:`msatk.stats.registry_extended`.
"""

from __future__ import annotations

from typing import Any

NUC = ("dna", "rna", "codon")
OBS = "Observed = canonical residue; gaps, missing, ambiguity codes and padding are unobserved."
PAIR_COST = "O(N^2 L) time, O(N^2) memory; skipped above max_exact_sequences."

Def = tuple[str, str, dict[str, Any]]

_Q = "Quantiles q05, q25, q50, q75, q95."

CHAR_DEFS: dict[str, list[Def]] = {
    "completeness": [
        (
            "alignment_completeness",
            "Fraction of all cells observed (AliStat Ca).",
            {"formula": "observed cells / (N L)", "value_range": "0-1", "missing_policy": OBS},
        ),
        (
            "min_sequence_completeness",
            "Lowest per-sequence completeness (AliStat Cr).",
            {"value_range": "0-1"},
        ),
        (
            "min_site_completeness",
            "Lowest per-column completeness (AliStat Cc).",
            {"value_range": "0-1"},
        ),
        ("sequence_completeness_quantiles", f"Per-sequence completeness. {_Q}", {}),
        ("site_completeness_quantiles", f"Per-column completeness. {_Q}", {}),
        ("min_pair_completeness_threshold", "Threshold used for pairs_below_threshold.", {}),
        ("pair_completeness_method", "exact, or why pairwise overlap was not computed.", {}),
        (
            "mean_pair_completeness",
            "Mean AliStat C_ij over pairs.",
            {"formula": "C_ij = columns observed in both i and j / L", "cost": PAIR_COST},
        ),
        ("min_pair_completeness", "Smallest C_ij: the least-overlapping pair.", {}),
        ("pair_completeness_quantiles", f"C_ij over all pairs. {_Q}", {}),
        ("pairs_below_threshold", "Pairs with C_ij below the threshold.", {}),
        ("pairs_without_overlap", "Pairs sharing no observed column (distance undefined).", {}),
        ("fraction_pairs_below_threshold", "pairs_below_threshold / pairs.", {}),
        ("min_shared_sites", "Fewest columns observed in both sequences of any pair.", {}),
    ],
    "sequence_completeness": [
        ("sequence_id", "Sequence identifier.", {}),
        ("sequence_completeness", "Observed columns / L (AliStat Cr).", {}),
        ("least_overlapping_partner", "Sequence with the smallest C_ij to this one.", {}),
        ("min_pair_completeness", "That smallest C_ij.", {}),
        ("pairs_below_threshold", "Partners with C_ij below the threshold.", {}),
    ],
    "gap_structure": [
        (
            "gap_openings",
            "Internal gap runs (a run starts after a non-gap cell).",
            {
                "missing_policy": "Terminal gaps (before the first / after the last non-gap cell) "
                "are counted separately, not as openings."
            },
        ),
        (
            "gap_extensions",
            "Internal gap cells beyond the first of each run.",
            {"formula": "internal_gap_cells - gap_openings"},
        ),
        ("internal_gap_cells", "Gap cells between a sequence's first and last non-gap cell.", {}),
        ("terminal_gap_cells", "Gap cells before the first / after the last non-gap cell.", {}),
        ("terminal_gap_fraction", "terminal_gap_cells / all gap cells.", {"undefined": "no gaps"}),
        ("mean_gap_run_length", "Mean internal gap-run length.", {"undefined": "no internal gaps"}),
        ("median_gap_run_length", "Median internal gap-run length.", {}),
        ("max_gap_run_length", "Longest internal gap run.", {}),
        ("single_site_gap_runs", "Internal gap runs of length 1.", {}),
        ("longest_gap_free_block", "Longest run of columns containing no gap.", {}),
        ("gap_free_blocks", "Number of maximal gap-free column runs.", {}),
        ("longest_complete_block", "Longest run of columns observed in every sequence.", {}),
        ("complete_blocks", "Number of maximal complete column runs.", {}),
        (
            "distinct_column_gap_patterns",
            "Distinct gap-presence patterns among gapped columns.",
            {},
        ),
        ("distinct_sequence_gap_patterns", "Distinct gap-position sets among sequences.", {}),
        (
            "mean_pairwise_gap_concordance",
            "Mean Jaccard similarity of gap positions over pairs.",
            {
                "formula": "|G_i ∩ G_j| / |G_i ∪ G_j| (1 when neither has gaps)",
                "value_range": "0-1",
            },
        ),
        (
            "fragmentation_index",
            "Mean over sequences of 1 - longest observed segment / observed "
            "cells; 0 when every sequence is one contiguous block.",
            {"value_range": "0-1"},
        ),
    ],
    "gap_structure_by_sequence": [
        ("sequence_id", "Sequence identifier.", {}),
        ("internal_gap_openings", "Internal gap runs.", {}),
        ("internal_gap_extensions", "Internal gap cells beyond each run's first.", {}),
        ("internal_gap_cells", "Internal gap cells.", {}),
        ("terminal_gap_cells", "Terminal gap cells.", {}),
        ("longest_internal_gap", "Longest internal gap run.", {}),
        ("observed_segments", "Maximal runs of observed cells.", {}),
        ("longest_observed_segment", "Longest such run.", {}),
        ("fragmentation", "1 - longest_observed_segment / observed cells.", {}),
        (
            "gap_pattern_distance",
            "1 - mean Jaccard similarity of this sequence's gap positions "
            "to every other sequence's.",
            {},
        ),
        (
            "gap_pattern_robust_z",
            "Robust z of gap_pattern_distance (median/MAD); large positive "
            "values mark unusual gap placement.",
            {"undefined": "MAD = 0"},
        ),
    ],
    "gap_run_length_distribution": [
        ("run_length", "Internal gap-run length.", {}),
        ("runs", "Runs of that length.", {}),
    ],
    "contiguous_blocks": [
        ("block_type", "gap_free (no gap in the column) or complete (all observed).", {}),
        ("start", "First column (1-based).", {}),
        ("end", "Last column (1-based, inclusive).", {}),
        ("length", "Columns in the block.", {}),
    ],
    "redundancy": [
        ("identity_threshold", "Identity at or above which two sequences are neighbours.", {}),
        (
            "min_overlap",
            "Neighbours' comparable sites must cover at least this fraction of "
            "the longer sequence's residues.",
            {},
        ),
        (
            "effective_sequences",
            "Effective number of sequences N_eff.",
            {
                "formula": "sum_i 1/n_i, n_i = neighbours of i (including i)",
                "value_range": "1-N",
                "cost": PAIR_COST,
                "assumptions": "Identity over sites where both sequences have a residue.",
            },
        ),
        ("effective_sequence_fraction", "N_eff / N.", {"value_range": "0-1"}),
        ("effective_sequences_per_site", "N_eff / L.", {}),
        ("mean_neighbours", "Mean n_i.", {}),
        ("max_neighbours", "Largest n_i.", {}),
        ("clusters", "Connected components of the neighbour graph (single linkage).", {}),
        ("largest_cluster", "Size of the largest cluster.", {}),
        ("singleton_clusters", "Clusters with one sequence.", {}),
        ("pairwise_identity_quantiles", f"Pairwise identity. {_Q}", {}),
    ],
    "sequence_weights": [
        ("sequence_id", "Sequence identifier.", {}),
        ("neighbours", "n_i.", {}),
        ("weight", "w_i = 1/n_i; weights sum to N_eff.", {}),
        ("cluster", "Cluster label (0-based, first-seen order).", {}),
        ("cluster_size", "Sequences in the cluster.", {}),
    ],
    "cluster_size_distribution": [
        ("cluster_size", "Cluster size.", {}),
        ("clusters", "Clusters of that size.", {}),
    ],
    "site_patterns": [
        ("distinct_site_patterns", "Distinct columns (all characters, gaps included).", {}),
        ("site_pattern_ratio", "distinct_site_patterns / L.", {"value_range": "0-1"}),
        ("repeated_pattern_sites", "L - distinct_site_patterns.", {}),
        ("max_pattern_multiplicity", "Most frequent column pattern's count.", {}),
        (
            "site_pattern_entropy",
            "Shannon entropy (bits) of the pattern frequency distribution.",
            {},
        ),
        (
            "mean_gini_simpson",
            "Mean over observed columns of 1 - sum p_c^2.",
            {"missing_policy": "Canonical residues only."},
        ),
        ("mean_effective_alphabet_size", "Mean 2^H over observed columns.", {}),
        (
            "mean_weighted_entropy",
            "Mean column entropy with sequence weights w_i.",
            {"undefined": "no exact pairwise identities"},
        ),
        ("weighted", "Whether sequence weights were available.", {}),
        ("sites_1_state", "Columns with one observed state.", {}),
        ("sites_2_states", "Columns with two observed states.", {}),
        ("sites_3_states", "Columns with three observed states.", {}),
        ("sites_4plus_states", "Columns with four or more observed states.", {}),
        ("singleton_minor_sites", "Variable columns whose non-majority states total 1.", {}),
        ("doubleton_minor_sites", "Variable columns whose non-majority states total 2.", {}),
        ("minor_count_3plus_sites", "Variable columns whose non-majority states total >= 3.", {}),
    ],
    "site_diversity": [
        ("site_index", "Column (1-based).", {}),
        ("simpson", "sum p_c^2.", {}),
        ("gini_simpson", "1 - sum p_c^2.", {}),
        ("effective_alphabet_size", "2^H.", {}),
        ("weighted_entropy", "Entropy (bits) with sequence weights.", {}),
        ("minor_state_count", "Observed residues not in the majority state.", {}),
        ("pattern_multiplicity", "Columns sharing this column's pattern.", {}),
        ("pattern_first_site", "First column with this pattern.", {}),
    ],
    "composition_heterogeneity": [
        (
            "rcv",
            "Relative composition variability (Phillips & Penny 2003).",
            {
                "formula": "sum_c sum_i |x_ic - mean_c| / (N L)",
                "missing_policy": "Canonical residues only.",
            },
        ),
        ("max_rcvt", "Largest per-sequence RCV contribution.", {}),
        ("sequences_tested", "Sequences with at least one residue.", {}),
        ("sequences_failing_composition_test", "Sequences with χ² p < alpha.", {}),
        ("fraction_failing_composition_test", "Failing / tested.", {}),
        ("alpha", "Significance level.", {}),
        (
            "mean_composition_jsd",
            "Mean Jensen–Shannon divergence (bits) of sequence composition "
            "from the alignment composition.",
            {"value_range": "0-1"},
        ),
        ("max_composition_jsd", "Largest such divergence.", {}),
    ],
    "sequence_composition_tests": [
        ("sequence_id", "Sequence identifier.", {}),
        ("rcvt", "Per-sequence RCV contribution: sum_c |x_ic - mean_c| / (N L).", {}),
        (
            "composition_chi2",
            "Pearson χ² of residue counts vs alignment frequencies.",
            {"assumptions": "Sites treated as independent; a screening test, as in IQ-TREE."},
        ),
        ("composition_df", "Degrees of freedom: observed states - 1.", {}),
        ("composition_p", "χ² p-value.", {}),
        ("composition_test_failed", "p < alpha.", {}),
        ("composition_jsd", "Jensen–Shannon divergence (bits) from the alignment composition.", {}),
    ],
    "symmetry_tests": [
        ("pair_method", "exact (all pairs) or sampled.", {}),
        ("pairs_tested", "Pairs tested.", {}),
        ("pairs_total", "N(N-1)/2.", {}),
        ("sample_seed", "Random seed when sampled.", {}),
        ("alpha", "Significance level.", {}),
        (
            "fraction_pairs_rejecting_bowker",
            "Pairs with Bowker p < alpha (rejecting symmetry: "
            "stationary, reversible, homogeneous evolution).",
            {"assumptions": "No multiple-testing correction; see Naser-Khdour et al. 2019."},
        ),
        (
            "fraction_pairs_rejecting_marginal",
            "Pairs with Stuart–Maxwell marginal symmetry p < alpha (rejecting stationarity).",
            {},
        ),
        ("fraction_pairs_rejecting_internal", "Pairs with internal symmetry p < alpha.", {}),
        ("min_bowker_p", "Smallest Bowker p-value.", {}),
        ("multiple_testing", "Multiple-testing note.", {}),
    ],
    "pair_symmetry_tests": [
        ("sequence_a", "First sequence.", {}),
        ("sequence_b", "Second sequence.", {}),
        ("sites", "Sites where both carry a canonical residue (divergence-matrix total).", {}),
        ("bowker", "Bowker statistic sum_{i<j} (n_ij - n_ji)^2 / (n_ij + n_ji).", {}),
        ("bowker_df", "Off-diagonal pairs with n_ij + n_ji > 0.", {}),
        ("bowker_p", "χ² p-value.", {}),
        ("marginal", "Stuart–Maxwell statistic d' V^-1 d on the first k-1 margins.", {}),
        ("marginal_df", "Rank of V.", {}),
        ("marginal_p", "χ² p-value.", {}),
        ("internal", "bowker - marginal.", {}),
        ("internal_df", "bowker_df - marginal_df.", {}),
        ("internal_p", "χ² p-value.", {"undefined": "internal_df <= 0"}),
    ],
    "windows": [
        ("window", "Window width (sites; a multiple of 3 for codon alignments).", {}),
        ("step", "Step between window starts.", {}),
        ("windows", "Number of windows.", {}),
        ("most_variable_window", "Window with the highest variable_fraction.", {}),
        ("gappiest_window", "Window with the highest gap_fraction.", {}),
        ("highest_entropy_window", "Window with the highest mean_entropy.", {}),
        ("windows_with_outlier_sequences", "Windows with at least one outlier sequence.", {}),
        ("outlier_z_threshold", "Robust-z threshold for window outliers.", {}),
        ("outlier_mad_floor", "Lower bound on the scaled MAD (identity units).", {}),
    ],
    "site_properties": [
        ("site_index", "Column (1-based).", {}),
        ("hydropathy_sd", "SD of Kyte–Doolittle hydropathy over residues in the column.", {}),
        ("volume_sd", "SD of residue volume (Zamyatnin 1972, Å^3).", {}),
        ("charge_variance", "Variance of net charge (K, R +1; D, E -1; H +0.1).", {}),
        ("polar_fraction", "Fraction of residues that are polar (STNQCYHKRDE).", {}),
        ("substitutions_vs_consensus", "Residues differing from the column consensus.", {}),
        ("conservative_substitutions", "Of those, BLOSUM62(residue, consensus) > 0.", {}),
        ("radical_substitutions", "Of those, BLOSUM62(residue, consensus) < 0.", {}),
        (
            "class_changes",
            "Of those, a change of residue class (hydrophobic, polar, positive, "
            "negative, special G/P/C).",
            {},
        ),
        ("radical_fraction", "radical / substitutions.", {}),
    ],
    "sequence_substitution_scores": [
        ("sequence_id", "Sequence identifier.", {}),
        (
            "mean_blosum62_to_consensus",
            "Mean BLOSUM62 score of residues vs the column consensus.",
            {},
        ),
        (
            "substitution_robust_z",
            "Robust z of that mean; strongly negative = unusually radical replacements.",
            {},
        ),
        ("radical_substitutions", "Radical substitutions in the sequence.", {}),
        ("conservative_substitutions", "Conservative substitutions in the sequence.", {}),
    ],
    "protein_properties": [
        ("substitutions_vs_consensus", "Residues differing from their column consensus.", {}),
        ("conservative_fraction", "BLOSUM62 > 0 share of those substitutions.", {}),
        ("radical_fraction", "BLOSUM62 < 0 share.", {}),
        ("neutral_fraction", "BLOSUM62 = 0 share.", {}),
        ("class_change_fraction", "Share changing residue class.", {}),
        ("median_hydropathy_sd", "Median column hydropathy SD.", {}),
        (
            "sequences_with_substitution_z_below_minus_3_5",
            "Sequences with substitution robust z < -3.5.",
            {},
        ),
    ],
    "codon_evolution": [
        ("genetic_code", "NCBI genetic code table.", {}),
        ("codon_sites", "Codon columns (L / 3).", {}),
        ("variable_codon_sites", "Codon columns with more than one complete, non-stop codon.", {}),
        (
            "synonymous_only_codon_sites",
            "Variable codon columns whose codons all encode one amino acid.",
            {},
        ),
        ("nonsynonymous_codon_sites", "Variable codon columns with more than one amino acid.", {}),
        (
            "effective_number_of_codons",
            "Wright's ENC from all complete codons.",
            {"value_range": "20-61", "undefined": "a synonymous family class is unobserved"},
        ),
        ("fourfold_sites", "Third positions whose consensus codon is four-fold degenerate.", {}),
        (
            "fourfold_nucleotide_diversity",
            "Mean pairwise difference per four-fold site among "
            "codons sharing the consensus first two bases.",
            {},
        ),
        ("dnds_method", "Estimator used for dN and dS.", {}),
        ("dnds_pair_method", "exact (all pairs) or sampled.", {}),
        ("dnds_pairs", "Pairs evaluated.", {}),
        ("dnds_pairs_defined", "Pairs with a defined dN/dS.", {}),
        (
            "median_pairwise_dn_ds",
            "Median pairwise dN/dS.",
            {"assumptions": "Counting estimate between pairs; not a phylogenetic selection test."},
        ),
        ("mean_pairwise_ds", "Mean pairwise dS.", {}),
        ("mean_pairwise_dn", "Mean pairwise dN.", {}),
        ("dnds_pairs_saturated", "Pairs with p >= 0.75 (Jukes–Cantor undefined).", {}),
    ],
    "pairwise_dnds": [
        ("sequence_a", "First sequence.", {}),
        ("sequence_b", "Second sequence.", {}),
        ("codons", "Codons compared (complete, non-stop in both).", {}),
        ("synonymous_sites", "S (Nei–Gojobori sites, averaged over the two sequences).", {}),
        ("nonsynonymous_sites", "N = 3 codons - S.", {}),
        ("synonymous_differences", "S_d.", {}),
        ("nonsynonymous_differences", "N_d.", {}),
        ("ps", "S_d / S.", {}),
        ("pn", "N_d / N.", {}),
        ("ds", "Jukes–Cantor corrected p_S.", {"formula": "-3/4 ln(1 - 4 p_S / 3)"}),
        ("dn", "Jukes–Cantor corrected p_N.", {}),
        ("dn_ds", "dN / dS.", {"undefined": "dS = 0, saturation, or < 10 codons"}),
        ("note", "Why a value is undefined.", {}),
    ],
    "coevolution": [
        (
            "candidate_sites",
            "Columns analysed: variable, occupancy >= min_occupancy, the max_sites most entropic.",
            {},
        ),
        ("max_sites", "Cap on candidate columns.", {"cost": "O(N (S k)^2) time."}),
        ("min_occupancy", "Minimum fraction of sequences observed in a candidate column.", {}),
        ("pseudocount", "lambda in p = (1 - lambda) f + lambda / q.", {}),
        ("weighting", "Sequence weighting used.", {}),
        ("effective_sequences", "Sum of weights.", {}),
        ("min_separation", "Minimum column distance for reported pairs.", {}),
        ("pairs", "Column pairs evaluated.", {}),
        ("note", "Why no pairs were evaluated.", {}),
        ("mean_mi", "Mean MI over evaluated pairs (bits).", {}),
        ("max_mi_apc", "Largest APC-corrected MI.", {}),
        ("top_pair", "Columns of the highest-ranked pair.", {}),
        ("pairs_reported", "Rows in coevolution_pairs.csv.", {}),
    ],
    "coevolution_pairs": [
        ("rank", "Rank by APC-corrected MI.", {}),
        ("site_i", "First column (1-based).", {}),
        ("site_j", "Second column (1-based).", {}),
        ("mi", "Mutual information (bits).", {"formula": "sum p(x,y) log2 p(x,y) / (p(x) p(y))"}),
        ("mi_apc", "MI - MI(i,.) MI(.,j) / MI(.,.) (Dunn et al. 2008).", {}),
        ("normalized_mi", "2 MI / (H_i + H_j).", {"value_range": "0-1"}),
        ("joint_entropy", "H(i, j) in bits.", {}),
        ("conditional_entropy_j_given_i", "H(j | i) = H(i, j) - H(i).", {}),
        ("effective_pair_weight", "Weighted sequences with residues in both columns.", {}),
    ],
    "comparison": [
        ("pair_method", "exact (all sequence pairs) or sampled.", {}),
        ("sequence_pairs", "Sequence pairs compared.", {}),
        ("reference_residue_pairs", "Aligned residue pairs implied by the reference.", {}),
        ("test_residue_pairs", "Aligned residue pairs implied by the test alignment.", {}),
        ("shared_residue_pairs", "Pairs implied by both.", {}),
        (
            "sp_score",
            "Sum-of-pairs score: shared / reference pairs (recall).",
            {"value_range": "0-1"},
        ),
        ("precision", "shared / test pairs.", {"value_range": "0-1"}),
        ("f1", "Harmonic mean of SP and precision.", {}),
        (
            "tc_score",
            "Total-column score: reference columns (>= 2 residues) reproduced exactly.",
            {"value_range": "0-1"},
        ),
        ("reference_columns", "Reference columns with at least two residues.", {}),
        ("mean_column_agreement", "Mean column agreement over test columns.", {}),
        ("columns_below_0_5_agreement", "Test columns with agreement < 0.5.", {}),
    ],
    "column_agreement": [
        ("site_index", "Test-alignment column (1-based).", {}),
        ("residue_pairs", "Residue pairs in the column (over compared sequence pairs).", {}),
        ("agreement", "Fraction of those pairs also aligned in the other alignment.", {}),
    ],
    "sequence_agreement": [
        ("sequence_id", "Sequence identifier.", {}),
        (
            "pairs_recovered",
            "Fraction of the reference's residue pairs involving this sequence "
            "that the test alignment recovers.",
            {},
        ),
    ],
    "ensemble": [
        ("alternatives", "Alternative alignments compared.", {}),
        ("mean_sp_vs_alternatives", "Mean SP of the test alignment against each alternative.", {}),
        ("mean_tc_vs_alternatives", "Mean TC against each alternative.", {}),
        ("mean_column_confidence", "Mean per-column confidence.", {}),
        ("columns_below_0_5_confidence", "Columns with confidence < 0.5.", {}),
        ("columns_fully_consistent", "Columns aligned identically in every alternative.", {}),
    ],
    "column_confidence": [
        ("site_index", "Test-alignment column (1-based).", {}),
        ("confidence", "Mean column agreement over the alternatives.", {}),
        ("min_agreement", "Lowest agreement with any alternative.", {}),
    ],
}

CHAR_DEFS["windows_table"] = [
    ("window_start", "First column (1-based).", {}),
    ("window_end", "Last column (inclusive).", {}),
    ("sites", "Columns in the window.", {}),
    ("occupancy", "Mean column occupancy.", {}),
    ("gap_fraction", "Mean column gap fraction.", {}),
    ("variable_fraction", "Fraction of variable columns.", {}),
    ("parsimony_informative_fraction", "Fraction of parsimony-informative columns.", {}),
    ("mean_entropy", "Mean column entropy (bits).", {}),
    (
        "mean_consensus_identity",
        "Mean identity to the consensus in the window, over sequences "
        "with at least half the window observed.",
        {},
    ),
    ("outlier_sequences", "Sequences with window consensus-identity robust z < threshold.", {}),
    ("gc_content", "GC fraction (nucleotides).", {"applies_to": NUC}),
]

# Table name -> registry level (table names that collide with a summary level).
TABLE_LEVELS = {"windows": "windows_table"}

CHAR_TITLES = {
    "completeness": "Completeness (summary)",
    "sequence_completeness": "Completeness (per sequence)",
    "gap_structure": "Gap structure (summary)",
    "gap_structure_by_sequence": "Gap structure (per sequence)",
    "gap_run_length_distribution": "Gap-run length distribution",
    "contiguous_blocks": "Contiguous blocks",
    "redundancy": "Redundancy and N_eff (summary)",
    "sequence_weights": "Sequence weights and clusters",
    "cluster_size_distribution": "Cluster size distribution",
    "site_patterns": "Site diversity and patterns (summary)",
    "site_diversity": "Site diversity (per site)",
    "composition_heterogeneity": "Composition heterogeneity (summary)",
    "sequence_composition_tests": "Composition tests (per sequence)",
    "symmetry_tests": "Symmetry tests (summary)",
    "pair_symmetry_tests": "Symmetry tests (per pair)",
    "windows": "Windowed statistics (summary)",
    "windows_table": "Windowed statistics (per window)",
    "protein_properties": "Protein properties (summary)",
    "site_properties": "Protein properties (per site)",
    "sequence_substitution_scores": "Protein substitution scores (per sequence)",
    "codon_evolution": "Codon evolution (summary)",
    "pairwise_dnds": "Pairwise dN/dS",
    "coevolution": "Coevolution (summary)",
    "coevolution_pairs": "Coevolution (top column pairs)",
    "comparison": "Alignment comparison (summary)",
    "column_agreement": "Alignment comparison (per column)",
    "sequence_agreement": "Alignment comparison (per sequence)",
    "ensemble": "Alignment ensemble (summary)",
    "column_confidence": "Alignment ensemble (per column)",
}

ALIGNMENT_DEFS: list[Def] = [
    (
        "characterization",
        "Nested summaries of the characterization analyses that ran (levels below).",
        {},
    ),
    (
        "characterization_profile",
        "Metrics profile used (core, comprehensive, evolutionary, "
        "coevolution, all, none or custom).",
        {},
    ),
    (
        "characterization_applicability",
        "Why each characterization analysis ran or was skipped.",
        {},
    ),
]
