"""Metric dictionary: one explicit definition for every statistic MSATK reports.

Each :class:`MetricDefinition` states the resolution (alignment / sequence / site /
pair), a description, the formula or rule, the denominator, how gaps / missing /
ambiguous characters are handled, the value range and the molecule types for
which it is reported. ``undefined`` documents when the metric is null.

The registry is exported as ``metric_dictionary.json`` with every profile run and
rendered to ``docs/metrics.md``; tests assert that every emitted field is defined.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

ALL = ("dna", "rna", "codon", "protein", "translated_cds", "mixed", "unknown")
NUC = ("dna", "rna", "codon")

GAPS_EXCLUDED = "Gaps, missing (?, N for nucleotides, X for proteins) and ambiguity codes excluded."
CANONICAL_ONLY = "Only canonical residues counted; gaps, missing and ambiguous cells excluded."
POPGEN = (
    "Computed on complete sites (every sequence carries a canonical base). Assumes homologous "
    "population samples; interpret cautiously for divergent multi-species alignments."
)


@dataclass(frozen=True)
class MetricDefinition:
    name: str
    level: str
    description: str
    formula: str = ""
    denominator: str = ""
    missing_policy: str = ""
    value_range: str = ""
    applies_to: tuple[str, ...] = ALL
    undefined: str = ""

    def as_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["applies_to"] = list(self.applies_to)
        return data


def _m(name: str, level: str, description: str, **kwargs: object) -> MetricDefinition:
    return MetricDefinition(name=name, level=level, description=description, **kwargs)  # type: ignore[arg-type]


_A, _S, _Q, _P = "alignment", "site", "sequence", "pair"

METRICS: tuple[MetricDefinition, ...] = (
    # ---------------------------------------------------------------- alignment
    _m("tool", _A, "Name of the generating tool (always 'MSATK')."),
    _m("number_of_sequences", _A, "Number of sequences (rows), N.", value_range=">= 1"),
    _m("alignment_length", _A, "Aligned width L (number of columns).", value_range=">= 1"),
    _m("total_cells", _A, "N x L matrix cells, including padding.", formula="N * L"),
    _m("rectangular_alignment", _A, "True when every input sequence has the same aligned length."),
    _m(
        "sequence_lengths_consistent",
        _A,
        "Alias of rectangular_alignment (kept for compatibility).",
    ),
    _m(
        "padded_cells",
        _A,
        "Cells added to shorter input rows to form a rectangle. Padding is recorded as missing "
        "data ('?'), never as alignment gaps.",
        value_range=">= 0",
    ),
    _m("detected_format", _A, "Input file format used for parsing."),
    _m("detected_file_format", _A, "Alias of detected_format."),
    _m("molecule_type", _A, "Molecule type used for analysis (detected or user-specified)."),
    _m("detected_molecule_type", _A, "Alias of molecule_type."),
    _m("length_divisible_by_three", _A, "True when L mod 3 == 0."),
    _m(
        "appears_codon_aware",
        _A,
        "Nucleotide alignment with L divisible by 3 and no frame warnings.",
        applies_to=NUC,
    ),
    _m("stop_codons_exist", _A, "Any in-frame stop codon (TAA/TAG/TGA) observed.", applies_to=NUC),
    _m(
        "stop_codon_count", _A, "Number of in-frame stop codons over all sequences.", applies_to=NUC
    ),
    _m("internal_stop_codons_exist", _A, "Any stop codon before the final codon.", applies_to=NUC),
    _m(
        "frameshift_warnings_exist",
        _A,
        "Any frame warning (L or ungapped length not divisible by 3, or partially gapped codons).",
        applies_to=NUC,
    ),
    _m("frameshift_warning_count", _A, "Number of frame warnings.", applies_to=NUC),
    _m("translated_cds_mode", _A, "True when analysed as translated CDS."),
    _m("duplicate_id_count", _A, "Number of sequence identifiers that occur more than once."),
    _m(
        "unique_sequence_count",
        _A,
        "Number of distinct aligned sequences (exact string identity after case normalisation).",
    ),
    _m("duplicate_sequence_count", _A, "Sequences identical to an earlier sequence."),
    _m("residue_cells", _A, "Cells holding a canonical residue.", missing_policy=CANONICAL_ONLY),
    _m("ambiguous_cells", _A, "Cells holding a partially informative ambiguity code."),
    _m("missing_cells", _A, "Cells holding missing data (excludes padding)."),
    _m("stop_cells", _A, "Cells holding '*' (protein alignments)."),
    _m("invalid_cells", _A, "Cells holding characters outside the molecule's alphabet."),
    _m("total_gaps", _A, "Number of gap cells ('-' or '.')."),
    _m(
        "gap_fraction",
        _A,
        "Fraction of all cells that are gaps.",
        formula="total_gaps / (N*L)",
        value_range="[0, 1]",
    ),
    _m(
        "occupancy",
        _A,
        "Fraction of all cells holding a canonical residue.",
        formula="residue_cells / (N*L)",
        value_range="[0, 1]",
    ),
    _m(
        "ambiguous_fraction",
        _A,
        "Fraction of non-padding cells holding ambiguity codes.",
        value_range="[0, 1]",
    ),
    _m("mean_sequence_length", _A, "Mean ungapped sequence length.", undefined="No sequences."),
    _m("median_sequence_length", _A, "Median ungapped sequence length."),
    _m("min_sequence_length", _A, "Minimum ungapped sequence length."),
    _m("max_sequence_length", _A, "Maximum ungapped sequence length."),
    _m("sequence_length_sd", _A, "Population standard deviation of ungapped lengths."),
    _m(
        "constant_sites",
        _A,
        "Columns with exactly one observed canonical state. All-gap/all-missing columns are NOT "
        "constant; they are counted as unobserved.",
        missing_policy=GAPS_EXCLUDED,
    ),
    _m("conserved_sites", _A, "Alias of constant_sites (kept for compatibility)."),
    _m(
        "variable_sites",
        _A,
        "Columns with >= 2 observed canonical states.",
        missing_policy=GAPS_EXCLUDED,
    ),
    _m("unobserved_sites", _A, "Columns with no canonical residue (statistics undefined)."),
    _m(
        "parsimony_informative_sites",
        _A,
        "Columns with >= 2 states each observed in >= 2 sequences.",
        missing_policy=GAPS_EXCLUDED,
    ),
    _m(
        "singleton_sites",
        _A,
        "Variable columns that are not parsimony-informative (MEGA definition).",
        missing_policy=GAPS_EXCLUDED,
    ),
    _m("gap_free_sites", _A, "Columns without any gap."),
    _m("complete_sites", _A, "Columns in which every sequence carries a canonical residue."),
    _m(
        "entropy_mean",
        _A,
        "Mean site Shannon entropy (bits) over observed columns.",
        missing_policy="Unobserved columns excluded.",
        undefined="No observed column.",
    ),
    _m("entropy_median", _A, "Median site entropy over observed columns."),
    _m(
        "normalized_entropy_mean",
        _A,
        "Mean of entropy / log2(K) over observed columns.",
        value_range="[0, 1]",
    ),
    _m(
        "mean_consensus_identity",
        _A,
        "Mean over sequences of identity to the majority consensus.",
        value_range="[0, 1]",
    ),
    _m(
        "pairwise_method",
        _A,
        "How pairwise identity was computed: exact (all pairs, full matrix), all_pairs (all pairs, "
        "no matrix), sampled (deterministic random sample with replacement) or none.",
    ),
    _m("pairwise_pairs_total", _A, "N(N-1)/2."),
    _m("pairwise_pairs_evaluated", _A, "Number of pairs actually compared."),
    _m(
        "pairwise_pairs_without_comparable_sites",
        _A,
        "Evaluated pairs sharing no comparable site (identity undefined, excluded from means).",
    ),
    _m("pairwise_sample_seed", _A, "Random seed used for sampled pairwise comparison."),
    _m(
        "mean_pairwise_identity",
        _A,
        "Mean pairwise identity over evaluated pairs with >= 1 comparable site.",
        formula="matches / comparable_sites per pair",
        missing_policy="A site is comparable only when both sequences carry a canonical residue.",
        value_range="[0, 1]",
        undefined="Pairwise comparison disabled, N < 2, or no pair has comparable sites.",
    ),
    _m("median_pairwise_identity", _A, "Median pairwise identity over evaluated pairs."),
    _m("min_pairwise_identity", _A, "Minimum pairwise identity over evaluated pairs."),
    _m("max_pairwise_identity", _A, "Maximum pairwise identity over evaluated pairs."),
    _m("mean_pairwise_distance", _A, "1 - mean_pairwise_identity (uncorrected p-distance)."),
    _m(
        "gc_content",
        _A,
        "Pooled G+C fraction over all canonical bases.",
        formula="(G + C) / (A + C + G + T)",
        missing_policy=CANONICAL_ONLY,
        applies_to=NUC,
        undefined="No canonical bases (always null for proteins).",
    ),
    _m(
        "gc_skew",
        _A,
        "(G - C) / (G + C) over all canonical bases.",
        value_range="[-1, 1]",
        applies_to=NUC,
    ),
    _m(
        "at_skew",
        _A,
        "(A - T) / (A + T) over all canonical bases.",
        value_range="[-1, 1]",
        applies_to=NUC,
    ),
    _m(
        "popgen_sites_used",
        _A,
        "Number of complete sites used for population-genetic statistics.",
        applies_to=NUC,
    ),
    _m(
        "segregating_sites",
        _A,
        "Number of polymorphic complete sites, S.",
        missing_policy=POPGEN,
        applies_to=NUC,
    ),
    _m(
        "nucleotide_diversity_pi",
        _A,
        "Nucleotide diversity: mean pairwise differences per complete site (Nei & Li 1979).",
        formula="sum_sites [C(N,2) - sum_c C(n_c,2)] / C(N,2) / sites_used",
        missing_policy=POPGEN,
        applies_to=NUC,
        undefined="N < 2 or no complete site.",
    ),
    _m(
        "theta_watterson",
        _A,
        "Watterson's theta per site.",
        formula="S / a1 / sites_used, a1 = sum_{i=1}^{N-1} 1/i",
        missing_policy=POPGEN,
        applies_to=NUC,
    ),
    _m(
        "tajimas_d",
        _A,
        "Tajima's D (Tajima 1989).",
        missing_policy=POPGEN,
        applies_to=NUC,
        undefined="N < 4 or S = 0.",
    ),
    _m(
        "haplotype_count",
        _A,
        "Distinct haplotypes over complete sites.",
        missing_policy=POPGEN,
        applies_to=NUC,
    ),
    _m(
        "haplotype_diversity",
        _A,
        "Nei's haplotype diversity over complete sites.",
        formula="N/(N-1) * (1 - sum_i p_i^2)",
        missing_policy=POPGEN,
        value_range="[0, 1]",
        applies_to=NUC,
    ),
    _m(
        "transitions_pairwise",
        _A,
        "Transition differences summed over all pairs at complete sites.",
        applies_to=NUC,
    ),
    _m(
        "transversions_pairwise",
        _A,
        "Transversion differences summed over all pairs at complete sites.",
        applies_to=NUC,
    ),
    _m(
        "ts_tv_ratio",
        _A,
        "Ratio of pooled pairwise transition to transversion differences (uncorrected).",
        applies_to=NUC,
        undefined="No transversions observed.",
    ),
    # ---------------------------------------------------------------- provenance
    _m("msatk_version", _A, "MSATK version that produced the output."),
    _m("output_schema_version", _A, "Version of the MSATK output schema."),
    _m("input", _A, "Input path as given."),
    _m("input_sha256", _A, "SHA-256 digest of the input file."),
    _m("generated_at", _A, "UTC timestamp of the run (ISO 8601)."),
    _m("python_version", _A, "Python interpreter version."),
    _m("platform", _A, "Operating system / platform string."),
    _m("coordinate_convention", _A, "Coordinate convention of site indices (1-based, inclusive)."),
    _m("interpretation", _A, "Plain-language interpretation of the profile."),
    _m("recommended_next_steps", _A, "Suggested follow-up actions."),
    # ---------------------------------------------------------------- site
    _m("site_index", _S, "1-based column index in the input alignment."),
    _m("site_class", _S, "constant, variable or unobserved (no canonical residue)."),
    _m(
        "observed_count",
        _S,
        "Sequences with a canonical residue in the column.",
        missing_policy=CANONICAL_ONLY,
    ),
    _m("occupancy", _S, "observed_count / N.", value_range="[0, 1]"),
    _m("gap_count", _S, "Gap cells in the column."),
    _m("gap_fraction", _S, "gap_count / N.", value_range="[0, 1]"),
    _m(
        "missing_fraction",
        _S,
        "Missing-data cells / N (gaps and ambiguity codes not included).",
        value_range="[0, 1]",
    ),
    _m("ambiguous_fraction", _S, "Ambiguity-code cells / N.", value_range="[0, 1]"),
    _m(
        "dominant_character",
        _S,
        "Most frequent canonical residue (ties broken by alphabet order).",
        undefined="Unobserved column.",
    ),
    _m(
        "dominant_character_frequency",
        _S,
        "Frequency of the dominant residue among observed residues.",
        formula="max_c n_c / observed_count",
        value_range="[1/K, 1]",
        undefined="Unobserved column.",
    ),
    _m(
        "conservation_score",
        _S,
        "Alias of dominant_character_frequency.",
        undefined="Unobserved column.",
    ),
    _m(
        "consensus_support",
        _S,
        "Fraction of ALL sequences carrying the dominant residue (penalises gaps).",
        formula="max_c n_c / N",
        value_range="[0, 1]",
    ),
    _m("observed_state_count", _S, "Number of distinct canonical residues in the column."),
    _m(
        "entropy",
        _S,
        "Shannon entropy in bits of the canonical residue distribution.",
        formula="-sum_c p_c log2 p_c",
        missing_policy=GAPS_EXCLUDED,
        value_range="[0, log2 K]",
        undefined="Unobserved column.",
    ),
    _m(
        "normalized_entropy",
        _S,
        "entropy / log2(K), K = canonical alphabet size.",
        value_range="[0, 1]",
    ),
    _m("effective_states", _S, "Effective number of states, 2^entropy.", value_range="[1, K]"),
    _m(
        "information_content",
        _S,
        "log2(K) - entropy (no small-sample correction).",
        value_range="[0, log2 K]",
    ),
    _m("constant", _S, "True when exactly one canonical state is observed."),
    _m("variable", _S, "True when >= 2 canonical states are observed."),
    _m("parsimony_informative", _S, "True when >= 2 states each occur in >= 2 sequences."),
    _m("singleton", _S, "True when variable but not parsimony-informative."),
    _m("gc_content", _S, "G+C fraction among canonical bases in the column.", applies_to=NUC),
    _m("codon_position", _S, "Codon position (1, 2, 3) of the column.", applies_to=("codon",)),
    # ---------------------------------------------------------------- sequence
    _m("row_index", _Q, "1-based row index in the input alignment."),
    _m("sequence_id", _Q, "Sequence identifier."),
    _m("original_length", _Q, "Length of the row as provided (including gaps), before padding."),
    _m("aligned_width", _Q, "Alignment width L."),
    _m("padded_cells", _Q, "Missing-data cells added to reach width L."),
    _m("ungapped_length", _Q, "Non-gap, non-padding characters in the row."),
    _m("residue_count", _Q, "Canonical residues in the row."),
    _m("span_start", _Q, "First non-gap column (1-based).", undefined="All-gap row."),
    _m("span_end", _Q, "Last non-gap column (1-based).", undefined="All-gap row."),
    _m("span_length", _Q, "span_end - span_start + 1 (0 for an all-gap row)."),
    _m("gap_count", _Q, "Gap cells in the row."),
    _m("gap_fraction", _Q, "gap_count / L.", value_range="[0, 1]"),
    _m("leading_gaps", _Q, "Gap cells before the first non-gap character."),
    _m("trailing_gaps", _Q, "Gap cells after the last non-gap character."),
    _m("internal_gap_count", _Q, "Gap cells between span_start and span_end."),
    _m("gap_run_count", _Q, "Number of contiguous gap runs (including terminal runs)."),
    _m("max_gap_run", _Q, "Length of the longest gap run."),
    _m("missing_count", _Q, "Missing-data cells in the row (padding excluded)."),
    _m("ambiguous_count", _Q, "Ambiguity-code cells in the row."),
    _m("ambiguous_fraction", _Q, "ambiguous_count / L."),
    _m("occupancy", _Q, "residue_count / L.", value_range="[0, 1]"),
    _m("unique_residue_count", _Q, "Distinct canonical residues observed in the row."),
    _m(
        "consensus_comparable_sites",
        _Q,
        "Columns where the row has a residue and the consensus is defined.",
    ),
    _m(
        "consensus_identity",
        _Q,
        "Fraction of comparable columns where the row matches the majority consensus.",
        value_range="[0, 1]",
        undefined="No comparable column.",
    ),
    _m(
        "composition_distance",
        _Q,
        "Total-variation distance between the row's residue composition and the pooled composition.",
        formula="0.5 * sum_c |p_c(row) - p_c(alignment)|",
        value_range="[0, 1]",
    ),
    _m(
        "duplicate_of",
        _Q,
        "Identifier of the first earlier row with an identical sequence, if any.",
    ),
    _m("gc_content", _Q, "G+C fraction among canonical bases.", applies_to=NUC),
    _m("gc_skew", _Q, "(G - C) / (G + C).", applies_to=NUC),
    _m("at_skew", _Q, "(A - T) / (A + T).", applies_to=NUC),
    _m(
        "mean_identity_to_others",
        _Q,
        "Mean identity to every other sequence (exact pairwise mode only).",
        undefined="Pairwise matrix not computed or no comparable partner.",
    ),
    _m("nearest_neighbor", _Q, "Most similar other sequence (exact pairwise mode only)."),
    _m("nearest_neighbor_identity", _Q, "Identity to the nearest neighbour."),
    _m("outlier_score", _Q, "1 - consensus_identity.", value_range="[0, 1]"),
    _m(
        "consensus_identity_robust_z",
        _Q,
        "Robust z-score of consensus_identity: (x - median) / (1.4826 * MAD).",
        undefined="Fewer than 3 defined values.",
    ),
    _m("composition_robust_z", _Q, "Robust z-score of composition_distance."),
    # ---------------------------------------------------------------- pair
    _m("sequence_a", _P, "First sequence of the pair."),
    _m("sequence_b", _P, "Second sequence of the pair."),
    _m("comparable_sites", _P, "Columns where both sequences carry a canonical residue."),
    _m("mismatches", _P, "Comparable columns with different residues."),
    _m(
        "identity",
        _P,
        "(comparable_sites - mismatches) / comparable_sites.",
        undefined="No comparable sites.",
    ),
    _m("p_distance", _P, "1 - identity (uncorrected).", undefined="No comparable sites."),
)

_INDEX = {(m.level, m.name): m for m in METRICS}


def get_metric(name: str, level: str = _A) -> MetricDefinition:
    return _INDEX[(level, name)]


def is_defined(name: str, level: str) -> bool:
    return (level, name) in _INDEX


def metric_dictionary() -> list[dict[str, object]]:
    return [m.as_dict() for m in METRICS]


def metrics_markdown() -> str:
    """Render the registry as Markdown (used for docs/metrics.md)."""

    lines = [
        "# Metric dictionary",
        "",
        "Generated from `msatk.stats.registry`. Every field MSATK writes is defined here.",
        "",
        "Character states: **residue** (canonical), **ambiguous** (IUPAC/B/Z/J), **missing** "
        "(`?`; `N` for nucleotides; `X` for proteins; padding), **gap** (`-`, `.`), **stop** (`*`).",
        "Undefined values are written as empty CSV cells / JSON `null`, never as zero.",
        "",
    ]
    for level in (_A, _Q, _S, _P):
        lines += [
            f"## {level.title()}-level",
            "",
            "| Field | Definition | Formula | Range | Undefined when | Applies to |",
            "|---|---|---|---|---|---|",
        ]
        for m in METRICS:
            if m.level != level:
                continue
            desc = m.description + (f" *{m.missing_policy}*" if m.missing_policy else "")
            applies = "all" if m.applies_to == ALL else ", ".join(m.applies_to)
            lines.append(
                f"| `{m.name}` | {desc} | {m.formula and '`' + m.formula + '`'} | {m.value_range} | {m.undefined} | {applies} |"
            )
        lines.append("")
    return "\n".join(lines)
