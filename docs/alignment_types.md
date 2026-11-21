# Alignment types

MSATK separates three things that are often conflated:

- the **alignment type**, meaning the data model: rows of one matrix, genomic blocks, mapped reads, graph paths, and so on;
- the **method** that produced it (global, local, semiglobal, profile-guided). MSATK records this when the format states it and never infers it from the final alignment;
- the **file format**. One FASTA file can hold a DNA, codon or protein alignment, and a SAM file holds mapped reads rather than a multiple sequence alignment.

Each input is loaded into the data model that matches it. Only analyses that are biologically and mathematically valid for that model are run. The full format matrix is in [Compatibility](compatibility.md), which is generated from the machine-readable registry.

```bash
msatk inspect input_file          # what the file is, its data model, and what can be computed
msatk formats                     # the compatibility matrix (Markdown or --format json)
msatk profile input_file          # dispatches to the right analysis for the data model
```

```python
from msatk import load_alignment, inspect_alignment

aln = load_alignment("genomes.maf.gz")
aln.kind          # 'block'
aln.summary()     # block, genome and rearrangement statistics
aln.tables()      # {'blocks': [...], 'genomes': [...]}
```

## Data models

| Kind | Model | Typical formats | What MSATK computes | What it deliberately does not compute |
|---|---|---|---|---|
| `sequence` | `SequenceAlignment` / `MSATK` | FASTA, PHYLIP, CLUSTAL, NEXUS, Stockholm, A2M/A3M, MSF, MEGA, SELEX, MASE, EMBOSS | Full statistics, QC, figures, report and ML encodings, plus the specialized analyses below | – |
| `block` | `BlockAlignment` | MAF, XMFA, AXT | Block-size distribution; per-genome aligned fraction, coverage and duplication; strand changes and collinearity breakpoints; per-block gap fraction, variability and identity | Concatenating blocks into one matrix |
| `read` | `ReadAlignment` | SAM, BAM, CRAM | Mapping rate, flags, MAPQ, clipping, indels, splicing, mismatch rate, base quality, insert size, depth and breadth of coverage | All-vs-all read identity; column entropy over a pileup |
| `mapping` | `MappingAlignment` | PAF, PSL, BLAST tabular, GAF | Hits, multi-hit queries, identity, alignment length, MAPQ or e-value, query and target coverage | Column-wise MSA statistics |
| `graph` | `GraphAlignment` | GFA 1, rGFA | Segments, links, components, tips, branching, paths and walks, and core/shell/private segments | Flattening into a rectangle |
| `profile_model` | `ProfileModel` | HMMER3, HHM | Model metadata; match-state emission entropy (HMMER3) | Sequence statistics |
| `structure` | planned | PDB, mmCIF | – | – |

## Specialized analyses for sequence alignments

Each analysis runs only when it applies. The reason it ran or was skipped is reported in `summary.json → analysis_applicability`, in `msatk inspect`, and in the report's *Specialized analyses* section.

| Analysis | Runs when | Output |
|---|---|---|
| Pairwise alignment | Exactly two sequences | Identity, similarity (proteins), mismatches, gap openings and extensions, terminal versus internal gaps, coverage, and a recomputed score. The score is compared with the score the file reports (EMBOSS). |
| Substitution scores | Protein | Mean BLOSUM62 score per column (`site_substitution_scores.csv`) and its summary |
| RNA secondary structure | Stockholm `#=GC SS_cons` on a nucleotide alignment | For each base pair: Watson–Crick, GU, non-canonical or incomplete; compensatory and consistent changes (`rna_base_pairs.csv`) |
| Profile states | A2M/A3M lowercase inserts, or Stockholm `#=GC RF` | Match-state coverage, deletions and inserted residues (`profile_states.csv`) |
| Partitions | NEXUS `charset` blocks or `--partitions` (RAxML/IQ-TREE style) | Occupancy, variability and entropy for each partition, plus a χ² test of compositional homogeneity (`partitions.csv`) |
| Genetic codes | Codon or nucleotide alignment | All NCBI translation tables, selected with `--genetic-code N`. They affect stop codons, frame QC, codon usage, RSCU synonymous families and translation. |

## Format notes

- **A3M and A2M are not ordinary FASTA.** Lowercase residues are insertions relative to match states. They are removed from the analysed matrix and counted separately. A2M writes `.` for insert gaps; A3M omits them.
- **MAF, XMFA and AXT store blocks.** Coordinates are converted to 0-based, half-open, forward-strand intervals. MAF minus-strand starts refer to the reverse complement and are converted. Reading a multi-block MAF with `read_alignment()` raises an error instead of concatenating the blocks.
- **SAM, BAM and CRAM** are summarised as mapped reads. Depth counts `M`, `=` and `X` bases of primary, non-duplicate reads, as `samtools depth` does by default. Reading BAM and CRAM needs `pysam`. CRAM may need `--reference`.
- **MEGA** `.` means "identical to the first sequence" and is expanded. In **MSF**, `.` and `~` are gaps.
- **Graphs** are analysed as graphs. GFA 2, GBZ and ODGI are recognised but must be converted to GFA 1 first.
- **Gzip:** every text format can be read from `.gz` files. Compression is detected from the file's content.
- Unaligned FASTA/FASTQ, VCF, BED and GFF/GTF are sequences, reads, variants or annotations. MSATK does not treat them as alignments.

## Support levels

Every format in the registry has one of these levels:

| Level | Meaning |
|---|---|
| stable | Tested and within the scope of the publication |
| beta | Implemented and tested on fixtures |
| experimental | Implemented, but there is no fixture in CI yet (CRAM) |
| planned | `msatk inspect` recognises it and explains what to do, but it is not readable yet |

A test checks that every stable or beta format has a fixture that loads into the declared data model. Another test checks that `docs/compatibility.md` matches the registry, so the documentation cannot claim support that has not been validated.
