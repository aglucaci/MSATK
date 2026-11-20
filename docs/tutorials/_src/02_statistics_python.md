# 2. Alignment statistics in Python

**You will:** compute site, sequence and pairwise statistics from Python, see how MSATK
treats gaps and undefined values, find outlier sequences and look up metric definitions.

**Data:** the Ras-family seed alignment (Pfam PF00071). Tables are pandas DataFrames
when pandas is installed and lists of dicts otherwise.

## Load and summarise

```python
from msatk import MSATK, datasets

msa = MSATK(datasets.path("ras_pf00071"), pairwise="exact")
s = msa.summary()
print(s["molecule_type"], s["number_of_sequences"], "x", s["alignment_length"])
print("occupancy           ", round(s["occupancy"], 3))
print("constant / variable ", s["constant_sites"], "/", s["variable_sites"])
print("parsimony-informative", s["parsimony_informative_sites"])
print("mean entropy (bits) ", round(s["entropy_mean"], 3))
```

Statistics are computed once from a vectorised encoding of the alignment and then
cached, so calling `summary()`, `sites()` and `sequences()` again costs nothing.

## Per-site statistics

```python
sites = msa.sites()
print(sites[["site_index", "site_class", "occupancy", "dominant_character",
             "entropy", "information_content"]].head(8).to_string(index=False))
```

`site_class` is one of `constant`, `variable` or `unobserved`. Information content is
measured over the residues that are present, so read it together with `occupancy`. Here
are the most conserved columns that at least 90% of sequences fill:

```python
well_filled = sites[sites["occupancy"] >= 0.9]
top = well_filled.sort_values("information_content", ascending=False).head(8)
print(top[["site_index", "dominant_character", "occupancy",
           "information_content"]].to_string(index=False))
```

The P-loop glycines and lysine (GxxxxGK[ST]) and the DTAG motif appear at the top, as
expected for a GTPase family.

## Per-sequence statistics and outliers

```python
seqs = msa.sequences()
cols = ["sequence_id", "ungapped_length", "gap_fraction", "consensus_identity",
        "consensus_identity_robust_z", "composition_robust_z", "nearest_neighbor"]
print(seqs.sort_values("consensus_identity_robust_z")[cols].head(5).to_string(index=False))
```

The robust z-scores are calculated as (x − median) / (1.4826 × MAD). They measure how far a
sequence's identity to the consensus, or its composition, is from the rest of the
alignment. Because they use the median and MAD, the outliers themselves barely move the
reference point. QC uses them to flag `DIVERGENT_SEQUENCES` and `COMPOSITION_OUTLIERS`.
They need only one pass over the alignment, so they also work for alignments too large
for all-vs-all comparison.

## Pairwise identity

```python
pairs = msa.pairwise()
print(len(pairs), "pairs")
print(pairs[["sequence_a", "sequence_b", "comparable_sites", "identity"]]
      .sort_values("identity").head(3).to_string(index=False))
```

Identity is computed over *comparable* sites, meaning columns where both sequences
have a residue. If two sequences share no such column, their identity is `None` rather
than 0, and the pair is left out of the mean.

```python
from msatk.models import Alignment, SequenceRecord
from msatk.stats import compute_statistics

toy = Alignment((SequenceRecord("a", "ACGT----"), SequenceRecord("b", "----ACGT"),
                 SequenceRecord("c", "ACGTACGT")))
bundle = compute_statistics(toy, "dna", pairwise="exact")
for row in bundle.pairwise.long_table():
    print(row["sequence_a"], row["sequence_b"], row["comparable_sites"], row["identity"])
print("pairs without comparable sites:", bundle.summary["pairwise_pairs_without_comparable_sites"])
print("mean identity (defined pairs only):", bundle.summary["mean_pairwise_identity"])
```

## Large alignments

Comparing every pair grows with N². The default `pairwise="auto"` compares all pairs up
to 1,000 sequences. Above that it compares a seeded random sample of pairs, and
`pairwise_method` records which was used. Use `pairwise="none"` to skip the comparison
entirely.

```python
fast = MSATK(datasets.path("ras_pf00071"), pairwise="none").summary()
print(fast["pairwise_method"], fast["mean_pairwise_identity"])
```

## What does a field mean?

Every output field is defined in the metric dictionary: what it measures, its formula,
how missing data is handled, its range, and when it is undefined.

```python
from msatk.stats.registry import get_metric

for name, level in [("information_content", "site"), ("outlier_score", "sequence"),
                    ("tajimas_d", "alignment")]:
    m = get_metric(name, level)
    print(f"{name}: {m.description}")
    if m.undefined:
        print(f"    undefined when: {m.undefined}")
```

`msatk metrics` prints the whole dictionary, which is also in the
[Metric dictionary](../metrics.md) page.
