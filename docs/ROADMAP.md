# MSATK Roadmap

Target: **MSATK v1.0 — one command, one complete alignment characterization,
multiple interoperable output formats, no hidden scientific assumptions**, and an
accompanying Application Note.

Status key: ✅ done in 0.2 · 🟡 partial · ⬜ not started

## P0 — publication-critical

| Item | Status |
|---|---|
| Correct statistical edge cases; unified gap/missing/ambiguity semantics | ✅ |
| Pairwise calculations and embeddings optional and size-aware | ✅ |
| Robust input validation and core metrics with a metric dictionary | ✅ |
| Reliable offline HTML reports (directory and single-file) | ✅ |
| Configurable alignment PNG/SVG/PDF figures and lossless ML tensors | ✅ |
| Stable CLI/API, versioned outputs and documentation | 🟡 schema 2.0; freeze at 1.0 |
| Validation on real and synthetic alignments | 🟡 synthetic + scikit-allel cross-check; real datasets pending |

## P1 — strong differentiators

| Item | Status |
|---|---|
| Interactive, searchable alignment viewer | ✅ |
| Offline sortable sequence/site tables | ✅ |
| Alignment-to-image windowing with reversible coordinate maps | ✅ |
| Parquet table output | ⬜ |
| PDF/TIFF figure export | ✅ |
| Performance benchmarks and memory-efficient processing | 🟡 harness + results; streaming for >1e9 cells pending |
| PyTorch/TensorFlow dataset adapters | ✅ (NumPy dataset + tensor adapters) |
| Mutual information / coevolution, sequence weighting, effective sequence count | ✅ MI + APC, N_eff; DCA not planned |
| BLOSUM-based protein similarity and conservation scores | ✅ BLOSUM62 column scores, pairwise similarity |
| Alternative genetic codes; synonymous/non-synonymous variation | ✅ all NCBI codes; Nei–Gojobori dN/dS, ENC |
| Completeness/overlap (AliStat), gap structure, composition heterogeneity, symmetry tests | ✅ |
| Alignment comparison (SP, TC, ensemble column confidence) | ✅ `msatk compare` |
| Tree-based diagnostics (treeness/RCV, concordance factors, homoplasy, saturation) | ⬜ needs a tree input |
| CAI (reference codon set), k-mer diversity | ⬜ |
| Comparison against AMAS, AliStat, CIAlign, trimAl on shared metrics | 🟡 AMAS done (paper analysis); AliStat pending |

## Alignment families (capability registry)

| Family | Priority | Status |
|---|---|---|
| Nucleotide, protein and codon MSAs | P0 | ✅ |
| Pairwise alignments (EMBOSS; PAF/PSL/BLAST hits) | P0 | ✅ beta |
| Specialized MSA formats (A2M, A3M, MSF, MEGA, SELEX, MASE) | P0 | ✅ beta |
| Gzip input for all text formats | P0 | ✅ |
| Machine-readable capability registry, `msatk inspect`, `msatk formats` | P0 | ✅ |
| Reference-based read alignments (SAM/BAM/CRAM QC) | P1 | 🟡 SAM tested; BAM needs pysam in CI; CRAM experimental |
| Whole-genome block alignments (MAF, XMFA, AXT) | P1 | ✅ beta; chain and HAL planned |
| RNA structure-aware alignments (SS_cons pairs, covariation) | P1 | ✅ beta |
| Profile alignments (A2M/A3M states, HMMER3, HHM) | P1 | ✅ beta; profile scoring planned |
| Multi-locus / partitioned alignments | P1 | ✅ beta |
| Graph and pangenome alignments (GFA 1, rGFA, GAF) | P2 | ✅ beta; GFA 2, GBZ, ODGI planned |
| Protein structure-aware alignments (PDB, mmCIF) | P2 | ⬜ recognised by `inspect` only |
| BLAST XML/JSON, SNP/genotype and morphological matrices | P2 | ⬜ |

## P2 — post-publication

- Visual alignment editing with undo and coordinate provenance
- Alignment-uncertainty / ensemble agreement statistics
- Additional biological annotations and a plugin/metric registry for third parties
- Large-scale tiled visualization and GPU acceleration

## Application Note checklist

- [x] Metric definitions and statistical conventions documented
- [x] Hand-calculable unit tests for edge cases
- [x] Independent validation of population-genetic statistics (scikit-allel)
- [x] Runtime/memory benchmark harness
- [ ] Shared-metric agreement with AMAS / AliStat on public alignments
- [ ] Use case 1: viral coding-sequence alignment
- [ ] Use case 2: protein-family alignment with ML tensor export
- [ ] Four-panel composite figure from real MSATK outputs
- [ ] Archive release (Zenodo DOI), PyPI and Bioconda 1.0
