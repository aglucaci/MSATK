# Tutorials

Hands-on walkthroughs of what MSATK can do. They are separate from the Application Note.
Every tutorial uses data bundled with the package (`msatk datasets`), so it runs offline
as written. The outputs shown were produced by running the code, and the docs build
re-runs them to keep them current.

| # | Tutorial | What you'll do | Data |
|---|---|---|---|
| 1 | [Profile an alignment in one command](01_quickstart.md) | `inspect`, `profile`, summary, QC findings, interactive report | Ras family (Pfam PF00071) |
| 2 | [Alignment statistics in Python](02_statistics_python.md) | Site, sequence and pairwise tables; undefined values; outliers; the metric dictionary | Ras family |
| 3 | [QC of a viral coding alignment](03_viral_codon_qc.md) | Reading frame, codon usage, codon-position entropy, QC-driven filtering, genetic codes | 400 influenza H3N2 HA sequences |
| 4 | [Protein families and publication figures](04_protein_family_figures.md) | BLOSUM62 and physicochemical conservation, alignment figures, logos, themes | Ras family |
| 5 | [Machine-learning representations](05_ml_encodings.md) | Lossless encodings, ambiguity policies, windows with coordinates, PyTorch | Ras family |
| 6 | [Alignments with extra structure](06_annotated_alignments.md) | RNA base-pair covariation, A3M profile states, partitions, pairwise alignments | Rfam tRNA (RF00005), toy files |
| 7 | [Beyond multiple sequence alignments](07_beyond_msas.md) | Genome blocks, mapped reads, PAF hits, pangenome graphs, profile HMMs | toy files |
| 8 | [Characterizing an alignment](08_characterization.md) | N_eff, pairwise overlap, windows, composition and symmetry tests, dN/dS, coevolution, alignment comparison | H3N2 HA, Ras, tRNA |

## Setup

```bash
pip install "msatk[plots,dataframes]"     # tutorial 5 also uses: pip install torch
```

## The bundled datasets

```bash
msatk datasets
```

Real datasets are public alignments, redistributed unchanged:
- Pfam (CC0);
- Rfam (CC0);
- NCBI GenBank (no restrictions on use).

The `toy_*` files are small synthetic examples, designed so that every number in the
tutorials can be checked by hand.
