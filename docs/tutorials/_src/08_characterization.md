# 8. Characterizing an alignment: reliability, information and evolution

**You will:** measure how much independent information an alignment holds, check
whether the overlap between sequences supports pairwise comparison, find local
problems with windows, test the compositional assumptions of substitution models,
summarise codon-level variation, detect coevolving columns, and score an alignment
against alternatives.

**Data:** the influenza H3N2 HA alignment (400 coding sequences), the Ras family and the
Rfam tRNA seed alignment.

## Profiles

`msatk profile` runs the *core* characterization by default: completeness, gap
structure and the effective number of sequences. Larger profiles add analyses:

```bash
msatk datasets h3n2_ha --out data
msatk profile data/h3n2_ha.fasta.gz --type codon --metrics evolutionary --plots none --no-alignment-figure --out h3n2
```

In Python, `characterize` runs a profile and `compute` runs only what the requested
metrics need:

```python
from msatk import MSATK, datasets

h3 = MSATK(datasets.path("h3n2_ha"), sequence_type="codon")
r = h3.characterize("evolutionary")
for name, why in r.applicability.items():
    print(f"{name:26} {why}")
```

## How many sequences are really there?

400 sequences of one influenza subtype are far from 400 independent observations.
The effective number of sequences counts each cluster of near-identical sequences
(here, ≥ 80 % identity) roughly once:

```python
red = r.summaries["redundancy"]
print({k: round(red[k], 4) for k in ["effective_sequences", "effective_sequence_fraction",
                                     "clusters", "largest_cluster"]})
print("identity quantiles:", {k: round(v, 3) for k, v in red["pairwise_identity_quantiles"].items()})
```

With N_eff ≈ 1.4 at 80 % identity, this alignment suits population-level analyses
(π, Tajima's D) and selection analyses within the subtype. It is not suitable for
methods that assume diverse, roughly independent sequences, such as coevolution or
protein-family profile building. The Ras family is the opposite case:

```python
ras = MSATK(datasets.path("ras_pf00071"))
print(ras.compute(["effective_sequence_fraction", "redundancy.effective_sequences"]))
```

## Do the sequences overlap enough?

Overall completeness can hide pairs of sequences that share few columns. AliStat-style
pairwise completeness `C_ij` is the fraction of columns observed in both sequences:

```python
comp = ras.characterize().summaries["completeness"]
for key in ["alignment_completeness", "min_pair_completeness", "min_shared_sites",
            "pairs_below_threshold"]:
    print(f"{key:24} {round(comp[key], 3)}")
gaps = ras.characterize().summaries["gap_structure"]
for key in ["gap_openings", "mean_gap_run_length", "max_gap_run_length",
            "terminal_gap_fraction", "longest_complete_block"]:
    print(f"{key:24} {round(gaps[key], 3)}")
```

## Local problems: windows

```python
windows = r.tables["windows"]
print(f"{len(windows)} windows of {r.summaries['windows']['window']} sites")
for w in windows[:6]:
    print(w["window_start"], w["window_end"], "variable", round(w["variable_fraction"], 2),
          "outlier sequences", w["outlier_sequences"])
```

Between 13 and 16 sequences are outliers in every window along the gene. Checked
against the bundled host labels (`h3n2_ha_hosts`), every one of them comes from a
non-human host (swine, canine, environmental or avian), the same group the QC in
[tutorial 3](03_viral_codon_qc.md) flags. The divergence runs the length of the gene,
which points to biology rather than to an alignment problem confined to one region.

## Are the model's compositional assumptions met?

```python
for name in ["composition_heterogeneity", "symmetry_tests"]:
    s = r.summaries[name]
    print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()
                 if k.startswith(("rcv", "fraction", "sequences_failing", "pairs_tested"))})
```

Few rejections at α = 0.05, and no more than chance would produce, mean a
stationary, homogeneous substitution model is reasonable here. The tRNA alignment is
different: it spans the tree of life, and its composition varies with the GC content
of each genome:

```python
trna = MSATK(datasets.path("trna_rf00005"))
t = trna.characterize("evolutionary").summaries
print("tRNA RCV", round(t["composition_heterogeneity"]["rcv"], 3),
      "| failing χ²:", t["composition_heterogeneity"]["sequences_failing_composition_test"],
      "| pairs rejecting marginal symmetry:",
      round(t["symmetry_tests"]["fraction_pairs_rejecting_marginal"], 3))
```

## Codon-level variation

```python
codon = r.summaries["codon_evolution"]
for key in ["variable_codon_sites", "synonymous_only_codon_sites", "nonsynonymous_codon_sites",
            "effective_number_of_codons", "median_pairwise_dn_ds", "dnds_pairs_saturated"]:
    print(f"{key:30} {round(codon[key], 3)}")
```

The median pairwise dN/dS is well below 1, as expected for a protein under purifying
selection overall. This is a counting estimate for pairs of sequences; to test for
selection at individual sites or branches, use codon models in HyPhy.

## Coevolving columns

Mutual information with the average product correction (APC) finds columns that
change together. In the tRNA alignment these should be the stems of the cloverleaf:

```python
from msatk.stats.rna_structure import base_pair_rows

pairs = {(p["pair_i"], p["pair_j"]) for p in base_pair_rows(trna.alignment, trna.annotations.ss_cons)}
co = trna.characterize("coevolution")
for p in co.tables["coevolution_pairs"][:8]:
    print(p["site_i"], p["site_j"], round(p["mi_apc"], 3),
          "base pair" if (p["site_i"], p["site_j"]) in pairs else "")
```

No structure was given to the MI calculation; it recovers the base pairs from the
sequences alone.

## How confident is the alignment?

Conservation does not tell you whether a column is aligned correctly; comparing
alignments does. Here the H3N2 alignment is compared with an alternative in which one
sequence has a codon-length gap at site 301, which shifts the rest of that sequence by
one codon. The other sequences get three trailing gap columns, so the alternative stays
rectangular and every sequence keeps exactly the same residues:

```python
from msatk.models import Alignment, SequenceRecord
from msatk.stats.consistency import compare_alignments

recs = h3.alignment.records
first = SequenceRecord(recs[0].id, recs[0].sequence[:300] + "---" + recs[0].sequence[300:])
others = tuple(SequenceRecord(r.id, r.sequence + "---") for r in recs[1:])
alt = Alignment((first,) + others)
summary, tables = compare_alignments(alt, h3.alignment)
print({k: round(summary[k], 4) for k in ["sp_score", "precision", "tc_score"]})
low = [c for c in tables["column_agreement"] if c["agreement"] is not None and c["agreement"] < 1]
print(f"{len(low)} columns with agreement < 1, sites {low[0]['site_index']}-{low[-1]['site_index']};",
      "lowest", round(min(c["agreement"] for c in low), 4))
```

The disagreement starts at site 304, the first column after the inserted gap, and runs
to the end. TC is the fraction of columns before the shift (300 of 1,701). Each affected column
loses only the pairs involving that one sequence, which is why agreement stays near 1:
with 400 sequences, one misplaced residue changes 399 of the 79,800 pairs in its column. With several alternative
alignments, `msatk compare --alternatives` turns this into a per-column confidence.

## What does a metric mean?

```bash
msatk metrics describe redundancy.effective_sequences
```
