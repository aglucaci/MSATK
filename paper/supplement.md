---
title: "Supplementary material — MSATK: comprehensive statistical profiling, visualization and machine-learning representations of multiple sequence alignments"
bibliography: references.bib
---

# S1 Data and reproducibility

| Dataset | Source | Size | Notes |
|---|---|---|---|
| PF00042 globin, PF00069 protein kinase, PF00071 Ras, PF00076 RRM | Pfam seed alignments via the InterPro API [@mistry2021pfam; @paysan2023interpro] | 37–73 sequences, 79–419 columns | Stockholm converted to FASTA; insert gaps (`.`) written as `-` |
| Influenza A(H3N2) haemagglutinin CDS | GenBank [@sayers2023genbank]; query `Influenza A virus[orgn] AND H3N2 AND hemagglutinin AND complete cds AND 2017[pdat]`, first 400 records, CDS via `efetch rettype=fasta_cds_na` | 400 × 1,701 nt | All CDS were 1,701 nt; used without re-alignment. Duplicates retained. Host labels parsed from record titles (`h3n2_ha_metadata.tsv`) |
| Synthetic DNA | `benchmarks/run_benchmarks.py::synthetic_fasta`, seed 1 | 1,000 and 5,000 × 10,000 | 5% random substitutions per site from a random ancestor; every tenth sequence carries a gap block of L/20 |

The scripts in `paper/analysis/` regenerate every number, table and figure
(`prepare_data.py`, `compare_tools.py`, `use_cases.py`, `make_figure1.py`);
machine-readable results are in `paper/analysis/results/`. Timings were measured on
an Apple M1 Pro (16 GB RAM, macOS) with Python 3.13, MSATK 0.2.0 and AMAS 1.0.

# S2 Agreement and runtime versus AMAS

**Table S1.** Full comparison (median of three runs; each run a separate process,
wall time and peak resident memory from `/usr/bin/time -l`). Missing (%) for AMAS is
"Missing_percent"; for MSATK it is (gap cells + missing-data cells) / total cells.
MSATK was run as `msatk stats --pairwise none --output-format json`; AMAS as
`AMAS.py summary`.

| Alignment | Type | N | L | Missing % (AMAS / MSATK) | Variable (AMAS / MSATK) | PI (AMAS / MSATK) | IUPAC-ambiguous cells | AMAS s / MB | MSATK s / MB |
|---|---|---|---|---|---|---|---|---|---|
| PF00042 globin | protein | 73 | 141 | 18.877 / 18.877 | 134 / 134 | 119 / 119 | 27 | 0.03 / 14 | 0.13 / 39 |
| PF00069 protein kinase | protein | 37 | 419 | 36.141 / 36.141 | 374 / 374 | 271 / 271 | 0 | 0.03 / 14 | 0.13 / 39 |
| PF00071 Ras | protein | 60 | 229 | 27.715 / 27.715 | 184 / 184 | 160 / 160 | 0 | 0.03 / 14 | 0.12 / 39 |
| PF00076 RRM | protein | 68 | 79 | 12.100 / 12.100 | 77 / 77 | 71 / 71 | 0 | 0.03 / 14 | 0.13 / 39 |
| Influenza H3N2 HA CDS | DNA | 400 | 1,701 | 0.000 / 0.000 | 789 / 789 | 666 / 666 | 22 | 0.17 / 28 | 0.19 / 62 |
| Synthetic | DNA | 1,000 | 10,000 | 0.500 / 0.500 | 10,000 / 10,000 | 10,000 / 10,000 | 0 | 2.53 / 210 | 1.04 / 297 |
| Synthetic | DNA | 5,000 | 10,000 | 0.500 / 0.500 | 10,000 / 10,000 | 10,000 / 10,000 | 0 | 12.34 / 950 | 4.20 / 1,153 |

MSATK uses more memory than AMAS (about 1.2× on the largest input) because it holds
the encoded N × L state arrays needed for per-site and per-sequence statistics.
AliStat [@wong2020alistat] could not be built in the benchmarking environment and
is therefore not included in Tables 2 and S1.

# S3 Feature coverage

**Table S2.** Feature coverage compiled from each tool's documentation (✓ supported;
◐ partial or via a separate step; — not provided). MSATK entries refer to v0.2.0.
*To be confirmed against the current releases of each tool before submission.*

| Capability | MSATK | AMAS | AliStat | CIAlign | trimAl | Jalview |
|---|---|---|---|---|---|---|
| Alignment summary statistics | ✓ | ✓ | ✓ | ◐ | ◐ | ◐ |
| Per-site statistics table | ✓ | — | ✓ | ◐ | ◐ | ◐ |
| Per-sequence statistics table | ✓ | — | ✓ | — | — | — |
| Pairwise identity / completeness | ✓ (exact or sampled) | — | ✓ | ✓ | ◐ | ✓ |
| Population-genetic summaries (π, θ<sub>W</sub>, Tajima's D) | ✓ | — | — | — | — | — |
| Codon-aware statistics (RSCU, frames, GC3) | ✓ | — | — | — | — | — |
| Explicit metric dictionary with undefined-value rules | ✓ | — | — | — | — | — |
| Structured QC findings | ✓ | — | — | ◐ | — | — |
| Alignment cleaning / trimming | — | — | ◐ | ✓ | ✓ | ◐ |
| Concatenation / format conversion | — | ✓ | — | — | ◐ | ✓ |
| Residue-level static figures and logos | ✓ | — | ◐ | ✓ | — | ✓ |
| Interactive viewer | ✓ (offline HTML) | — | — | — | — | ✓ (desktop) |
| ML tensors with manifests and windows | ✓ | — | — | — | — | — |
| Machine-readable provenance per run | ✓ | — | — | ◐ | — | — |

# S4 Use case 1: influenza H3N2 haemagglutinin

**Table S3.** QC findings (`msatk profile --type codon --preset comprehensive`).

| Code | Severity | Count |
|---|---|---|
| DIVERGENT_SEQUENCES | warning | 27 |
| COMPOSITION_OUTLIERS | warning | 21 |
| DUPLICATE_SEQUENCES | info | 134 |

**Table S4.** DIVERGENT_SEQUENCES and COMPOSITION_OUTLIERS by host label (from the
GenBank record title).

| Host label | Sequences | Divergent | Composition outlier |
|---|---|---|---|
| human | 380 | 11 | 9 |
| swine | 10 | 6 | 2 |
| canine | 4 | 4 | 4 |
| environment | 3 | 3 | 3 |
| avian (mallard, gadwall, blue-winged teal) | 3 | 3 | 3 |

Eighteen of the 21 composition outliers were also flagged as divergent. The four
swine sequences that were not flagged (MF280178.1, MF471676.1, MF471685.1,
MF471687.1) had identity to consensus of 0.986–0.989; their nearest neighbours were
a human strain (99.94%) or each other (99.7–99.8%).

**Table S5.** Population-genetic summaries before and after removing the 27
QC-flagged divergent sequences (complete sites; exact pairwise identity).

| Statistic | All 400 | After QC (373) |
|---|---|---|
| Segregating sites | 789 | 387 |
| Nucleotide diversity π | 0.0233 | 0.0119 |
| Watterson's θ (per site) | 0.0701 | 0.0347 |
| Tajima's D | −2.052 | −2.012 |
| Haplotype diversity | 0.992 | 0.991 |
| Ts/Tv | 2.85 | 3.45 |
| Mean pairwise identity | 0.976 | 0.988 |

Codon-position mean entropy (bits): position 1, 0.042; position 2, 0.038;
position 3, 0.142.

# S5 Use case 2: Ras family (PF00071)

**Table S6.** Guanine-nucleotide-binding motifs located by pattern search in
RAS2_YEAST and their per-site statistics. The G2 motif (a single conserved
threonine) was not searched.

| Motif | Pattern | Reference match | Alignment sites | Mean IC (bits) | Mean gap fraction |
|---|---|---|---|---|---|
| G1 (P-loop) | `G..G.GK[ST]` | GGGGVGKS | 6–13 | 3.48 | 0.00 |
| G3 | `D[TS]AG[QL]E` | DTAGQE | 69–74 | 4.25 | 0.00 |
| G4 | `[NT]K.D` | NKSD | 138–141 | 3.39 | 0.00 |
| G5 | `[ETS]A[KR]` | SAK | 209–211 | 3.86 | 0.00 |

Median information content of the other 141 columns with ≤ 50% gaps: 1.96 bits.
Median per-column standard deviation of Kyte–Doolittle hydropathy: motif columns
0.01, other occupied columns 1.79. A categorical encoding of the alignment decoded
to the input exactly. The windowed one-hot export (64 sites × 32 rows, stride 32)
produced 14 tiles of shape 32 × 64 × 22 (2 tiles with site padding and 7 with row
padding, all marked in the mask), written in 15 ms.
