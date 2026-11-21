# Statistical conventions

MSATK reports statistics at four resolutions — **alignment**, **sequence**, **site**
and **pair** — computed from one shared vectorised encoding so that every metric
treats characters identically. Every field is defined in the
[metric dictionary](metrics.md). Analyses beyond these base statistics (completeness,
N_eff, heterogeneity and symmetry tests, windows, coevolution, dN/dS, alignment
comparison) are described in [Alignment characterization](characterization.md).

## Character states

Each aligned cell is assigned one state, and the meaning depends on molecule type:

| State | Nucleotide (DNA/RNA/codon) | Protein |
|---|---|---|
| residue | `A C G T` (`U` ≡ `T`; RNA uses `A C G U`) | 20 standard amino acids (including **N** = asparagine) |
| ambiguous | `R Y S W K M B D H V` | `B Z J U O` |
| missing | `N`, `?` | `X`, `?` |
| gap | `-`, `.` | `-`, `.` |
| stop | — | `*` |

Earlier versions treated `N` and `X` as missing for every molecule type, which
silently discarded asparagine from protein statistics. Input is case-folded.

## Key rules

* **Unobserved columns.** A column with no canonical residue is *unobserved*. Its
  entropy and conservation are undefined (`null`) and it is not counted as a
  constant/conserved site.
* **Comparable sites.** Pairwise identity is `matches / comparable sites`, where a
  site is comparable only when both sequences carry a canonical residue. Pairs with
  no comparable site have undefined identity and are excluded from averages; their
  number is reported (`pairwise_pairs_without_comparable_sites`).
* **Unequal lengths.** Non-rectangular input is padded with missing data (`?`), not
  gaps, the number of padded cells is reported, and a QC warning is raised.
  `--strict` rejects such input.
* **Entropy** is Shannon entropy in bits over canonical residues only.
  `normalized_entropy = H / log2 K`.
* **Parsimony-informative** sites have ≥ 2 states each in ≥ 2 sequences; **singleton**
  sites are variable but not parsimony-informative.
* **GC content** at alignment level is pooled over all canonical bases.

## Pairwise cost policy

Exact all-vs-all identity is O(N²·L). MSATK computes it with blocked matrix
products over one-hot indicators. With `--pairwise auto` (default) the exact
matrix is used up to 1,000 sequences; above that, a deterministic random sample of
20,000 pairs (seeded, with replacement) estimates the identity distribution and the
method is recorded as `sampled`. `--pairwise none` skips comparison entirely.

## Population-genetic statistics

For nucleotide alignments MSATK reports segregating sites, nucleotide diversity π,
Watterson's θ, Tajima's D, haplotype count/diversity and pooled pairwise Ts/Tv.
These assume homologous population samples and are computed on **complete sites**
only (every sequence carries a canonical base); `popgen_sites_used` records how many.
π, θ_W and Tajima's D agree with scikit-allel to floating-point precision (see
`tests/unit/test_diversity.py`). Interpret them cautiously for divergent
multi-species alignments.

## Outlier scores

Sequence outliers use O(N·L) statistics that do not require a pairwise matrix:
identity to the majority consensus and total-variation distance between a
sequence's residue composition and the pooled composition, each converted to a
robust z-score (median/MAD).
