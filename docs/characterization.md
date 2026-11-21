# Alignment characterization

Beyond the base statistics (gaps, conservation, entropy, pairwise identity,
composition, popgen summaries), MSATK runs **characterization analyses**. Each one
answers a question about the alignment:

1. **Is the alignment reliable?** Completeness and pairwise overlap, gap structure,
   window outliers, and agreement with other alignments (`msatk compare`).
2. **How much information does it hold?** Effective number of sequences, site
   patterns, weighted diversity.
3. **What patterns are present?** Composition heterogeneity, protein property
   variation, codon-level selection summaries, coevolving columns.
4. **Is it suitable for downstream analysis?** Symmetry (stationarity) tests,
   overlap between sequence pairs, redundancy.

Each analysis runs only when it is valid for the input, and records why it ran or
was skipped (`characterization_applicability` in `summary.json`). Every field is
defined in the [metric dictionary](metrics.md), with its formula, missing-data
policy, assumptions and cost.

## Profiles

```bash
msatk profile alignment.fasta                          # core (default)
msatk profile alignment.fasta --metrics comprehensive
msatk profile alignment.fasta --metrics evolutionary
msatk profile alignment.fasta --metrics coevolution
msatk profile alignment.fasta --metrics all
```

| Analysis | core | comprehensive | evolutionary | coevolution | Runs when |
|---|:-:|:-:|:-:|:-:|---|
| `completeness` | ✓ | ✓ | ✓ | ✓ | always; pairwise overlap needs N ≤ `max_exact_sequences` |
| `gap_structure` | ✓ | ✓ | ✓ | ✓ | always |
| `redundancy` (N_eff) | ✓ | ✓ | ✓ | ✓ | N ≥ 2 and N ≤ `max_exact_sequences` |
| `site_patterns` |  | ✓ | ✓ | ✓ | always (weighted when N_eff is available) |
| `composition_heterogeneity` |  | ✓ | ✓ | ✓ | N ≥ 2 |
| `windows` |  | ✓ | ✓ | ✓ | L ≥ 20 |
| `protein_properties` |  | ✓ | ✓ | ✓ | protein alignments |
| `symmetry_tests` |  |  | ✓ |  | N ≥ 2 |
| `codon_evolution` |  |  | ✓ |  | codon alignments |
| `coevolution` |  |  |  | ✓ | N ≥ 20 |

`--preset comprehensive` uses the comprehensive profile; `basic` and `standard` use
core. Sequence divergence from the consensus (with robust z-scores), nucleotide
diversity, Watterson's θ, Tajima's D and haplotype diversity are base statistics
and are always reported.

## What each analysis computes

**Completeness** follows AliStat (Wong et al. 2020): a cell is observed when it
holds a canonical residue. MSATK reports the completeness of the whole alignment,
of each sequence and of each column, and **pairwise completeness** C_ij (columns
observed in both sequences / L). Two alignments with the same missing fraction can
differ greatly in C_ij, and C_ij determines whether a distance between two
sequences is supported by enough data. Tables: `sequence_completeness.csv`,
`pair_completeness_matrix.csv` (N ≤ 500).

**Gap structure** separates internal gaps (indel-like runs) from terminal gaps
(unsequenced ends). It reports openings, extensions, the run-length distribution,
the longest gap-free and complete blocks, gap-pattern diversity, gap-position
concordance between sequences, a fragmentation index, and a per-sequence
gap-pattern robust z-score that flags unusual gap placement.

**Redundancy**: N_eff = Σ 1/n_i, where n_i counts the sequences within 80 % identity
of sequence i whose comparable sites cover at least half of the longer sequence.
Reported with N_eff/N, N_eff/L, identity quantiles and single-linkage clusters at
the threshold. The weights w_i = 1/n_i are reused by the weighted site entropy and by
mutual information.

**Site patterns**: distinct column patterns and their multiplicity, pattern entropy,
Gini–Simpson diversity, effective alphabet size (2^H), sequence-weighted entropy,
and the distribution of sites by number of states and minor-state count.

**Composition heterogeneity**: RCV (Phillips & Penny 2003, as in PhyKIT) with
per-sequence contributions, a per-sequence χ² test against alignment frequencies
(as in IQ-TREE), and Jensen–Shannon divergence from the alignment composition.

**Symmetry tests**: for each pair of sequences (all pairs, or 20,000 sampled with a
recorded seed), Bowker's matched-pairs test, the Stuart–Maxwell marginal symmetry
test and the internal symmetry test (Naser-Khdour et al. 2019). Rejecting them
suggests the pair did not evolve under a stationary, reversible, homogeneous model.
No multiple-testing correction is applied; the fractions describe the uncorrected
tests.

**Windows**: occupancy, gaps, variability, entropy, GC and mean identity to the
consensus in sliding windows (default about L/20, a multiple of 3 for codons), plus
the number of sequences whose identity inside the window is an outlier (robust z <
−3.5).

**Protein properties**: per column, the SD of hydropathy and volume, the variance of
charge and the polar fraction. Substitutions relative to the consensus are classed
as conservative (BLOSUM62 > 0), radical (< 0) or neutral, and as changes of residue
class. A per-sequence BLOSUM62 robust z-score flags unusually radical replacements.

**Codon evolution**: pairwise dN, dS and dN/dS by Nei–Gojobori (1986) with
Jukes–Cantor correction under the selected genetic code. This is a counting estimate
for pairs of sequences, not a test of selection; use codon models in HyPhy or PAML
for that. Also reported: Wright's ENC, synonymous-only versus non-synonymous
variable codon sites, and nucleotide diversity at four-fold degenerate sites.

**Coevolution**: mutual information between columns with sequence weighting,
pseudocounts, the average product correction (APC) and normalised MI, computed on up
to 200 variable, well-occupied candidate columns and written to
`coevolution_pairs.csv`. This is a screening statistic; direct-coupling models are
out of scope.

## Choosing metrics in Python

```python
from msatk import MSATK

msa = MSATK("alignment.fasta")
r = msa.characterize("comprehensive")        # summaries, tables, applicability
r.summaries["redundancy"]["effective_sequences"]
r.tables["windows"]

msa.compute(["effective_sequence_fraction", "rcv", "sequence_weights",
             "redundancy.effective_sequences"])   # runs only what is needed
```

`compute` accepts metric names, `level.name` when a name exists at several levels,
analysis names (returning the analysis summary) and table names (returning rows).
`msatk metrics list` and `msatk metrics describe NAME` show what is available.

## Alignment confidence: `msatk compare`

High conservation does not mean a column is aligned correctly. Alignment accuracy
can only be estimated by comparing alignments, so MSATK keeps this separate from its
descriptive QC scores:

```bash
# Against a trusted reference (e.g. a BAliBASE benchmark)
msatk compare mafft.fasta --reference reference.fasta --out cmp/

# Against alternative alignments of the same sequences (other aligners/settings)
msatk compare mafft.fasta --alternatives muscle.fasta clustalo.fasta --out conf/
```

The reference comparison reports the **SP score** (fraction of reference residue
pairs recovered), precision, F1, the **TC score** (reference columns reproduced
exactly) and per-column and per-sequence agreement. With alternatives, the per-column
**confidence** is the mean agreement over the alternatives, in the spirit of GUIDANCE
column scores. Sequences are matched by ID, and their ungapped residues must be
identical across the alignments.

## Validation

The test suite checks:

* completeness, gap, N_eff, site-pattern, RCV and symmetry values on small alignments
  worked out by hand;
* dN and dS against Biopython's `NG86` (identical when MSATK is switched to
  Biopython's handling of stop codons; see the module documentation);
* that mutual information with APC recovers the tRNA cloverleaf: at least 9 of the
  top 10 pairs are base pairs of the Rfam RF00005 consensus structure;
* SP/TC on identical and shifted alignments;
* that every emitted field is defined in the metric dictionary.

## Not yet implemented

Statistics that need a tree or a model fit (treeness/RCV, site concordance factors,
homoplasy indices, saturation plots, relative rates), direct-coupling analysis,
CAI (which needs a reference codon set) and k-mer diversity are not computed. See
the [roadmap](ROADMAP.md).
