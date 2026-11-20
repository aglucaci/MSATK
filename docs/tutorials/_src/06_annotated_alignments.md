# 6. Alignments with extra structure: RNA, profiles, partitions and pairs

**You will:** use the information that some alignment files carry beyond the residues: a
consensus RNA secondary structure, match and insert states, partitions into genes, and the
method and score of a pairwise alignment. MSATK runs each of these analyses only when the
input supports it, and records why.

## RNA secondary structure and covariation

The Rfam tRNA seed alignment has a consensus cloverleaf structure in its
`#=GC SS_cons` line. Base pairs that are conserved *as pairs* while the individual bases
change (compensatory substitutions) are the classic evidence for a shared RNA structure.

```python
from msatk import load_alignment, datasets

trna = load_alignment(datasets.path("trna_rf00005"))
print(trna.alignment.n_sequences, "sequences x", trna.alignment.length, "columns;",
      "annotations:", trna.annotations.present())
print(trna.annotations.ss_cons)
```

```python
summary = trna.summary(pairwise="none")
for key, value in summary["specialized"]["rna_structure"].items():
    print(f"{key:40} {value}")
```

The table for each base pair shows which pairs covary most:

```python
pairs = trna.tables(pairwise="none")["rna_base_pairs"]
pairs.sort(key=lambda r: -r["compensatory_sequences"])
print(f"{'pair':>9} {'WC':>4} {'GU':>4} {'other':>5} {'compensatory':>12}  dominant")
for r in pairs[:6]:
    print(f"{r['pair_i']:>4}-{r['pair_j']:<4} {r['watson_crick']:>4} {r['wobble']:>4} "
          f"{r['noncanonical']:>5} {r['compensatory_sequences']:>12}  {r['dominant_pair']}")
```

The same alignment profiled from the command line writes these results to
`tables/rna_base_pairs.csv`, and the report gets a *Specialized analyses* section.

## Profile states (A2M/A3M)

In A2M and A3M files, lowercase letters are insertions relative to the profile's match
states. MSATK removes them from the matrix it analyses and counts them separately, so
insertions are never treated as aligned columns.

```python
from pathlib import Path

Path("family.a3m").write_text(">ref\nMKVLAAG\n>s1\nMKVaaLAAG\n>s2\nMK-LAAGkk\n>s3\nMRVLA-G\n")
a3m = load_alignment("family.a3m")
print(a3m.alignment.sequences)
prof = a3m.summary(pairwise="none")["specialized"]["profile_states"]
for key, value in prof.items():
    print(f"{key:40} {value}")
for row in a3m.tables(pairwise="none")["profile_states"]:
    print(row)
```

Stockholm files that mark match columns with a `#=GC RF` line are handled the same way.

## Partitioned (multi-locus) alignments

Concatenated matrices for phylogenomics are usually divided into genes or codon
positions. MSATK reads NEXUS `charset` blocks, or a RAxML/IQ-TREE partition file passed
with `--partitions`. It then reports each partition separately, including a χ² test of
whether base composition differs between sequences.

```python
nex = load_alignment(datasets.path("toy_partitioned"))
print([(p.name, p.sites) for p in nex.annotations.partitions])
for row in nex.tables(pairwise="none")["partitions"]:
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()
           if k in ("partition", "sites_in_alignment", "gap_fraction", "variable_sites",
                    "composition_chi2", "composition_chi2_p")})
```

```bash
msatk datasets toy_partitioned --out data
msatk profile data/toy_partitioned.nexus --type dna --out parts --plots none --no-alignment-figure
```

## Pairwise alignments

When an alignment has exactly two sequences, MSATK reports the statistics that pairwise
aligners report. If the file states the method and score, as EMBOSS output does, MSATK
keeps them and recomputes the score under the stated scheme, so the two can be compared.

```python
pair = load_alignment(datasets.path("toy_pairwise"))
print(pair.alignment.sequences, "| scope:", pair.scope, "| method:", pair.annotations.method)
stats = pair.summary()["specialized"]["pairwise_alignment"]
for key in ["identity", "identity_over_columns", "similarity", "gap_openings",
            "gap_extensions", "coverage_1", "coverage_2", "score_recomputed",
            "score_reported", "scoring_scheme"]:
    print(f"{key:24} {stats[key]}")
```

## Why did an analysis run or not?

```bash
msatk datasets trna_rf00005 --out data
msatk inspect data/trna_rf00005.sto.gz
```
