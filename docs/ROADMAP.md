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
| Mutual information / coevolution, sequence weighting, effective sequence count | ⬜ |
| BLOSUM-based protein similarity and conservation scores | ⬜ |
| Alternative genetic codes; synonymous/non-synonymous variation | ⬜ |
| Comparison against AMAS, AliStat, CIAlign, trimAl on shared metrics | ⬜ |

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
