# Benchmarks

Runtime and peak resident memory measured with `benchmarks/run_benchmarks.py`
(each configuration in a fresh process; synthetic DNA alignments derived from a
random ancestor with 5 % substitutions and gapped fragments in 10 % of rows).
Machine: Apple M1 Pro, 16 GB RAM, Python 3.13, NumPy 2.4 — MSATK 0.2.0.
Raw results: `benchmarks/results.csv`.

| Sequences × sites | Cells | Stats, no pairwise | Stats, `auto` pairwise | Stats, exact pairwise | `profile --preset basic` |
|---|---|---|---|---|---|
| 100 × 1,000 | 0.1 M | 0.11 s / 45 MB | 0.12 s / 46 MB (exact) | 0.12 s / 47 MB | 0.26 s / 52 MB |
| 100 × 10,000 | 1 M | 0.23 s / 91 MB | 0.26 s / 114 MB (exact) | 0.25 s / 109 MB | 0.45 s / 107 MB |
| 1,000 × 1,000 | 1 M | 0.22 s / 78 MB | 0.33 s / 196 MB (exact) | 0.31 s / 199 MB | 0.37 s / 91 MB |
| 1,000 × 10,000 | 10 M | 1.03 s / 322 MB | 1.22 s / 534 MB (exact) | 1.18 s / 536 MB | 1.22 s / 334 MB |
| 5,000 × 1,000 | 5 M | 0.73 s / 233 MB | 0.94 s / 250 MB (sampled) | 2.45 s / 1.34 GB | 1.37 s / 248 MB |
| 5,000 × 10,000 | 50 M | 4.37 s / 1.23 GB | 4.70 s / 1.21 GB (sampled) | 8.64 s / 2.01 GB | 5.25 s / 1.31 GB |

"Stats" computes every alignment-, sequence- and site-level statistic (including
population-genetic statistics). Times include Python start-up and file parsing.
Figures (matplotlib) are excluded because their cost depends on the formats requested.

Observations:

* Statistics scale linearly with the number of cells (~90 ns/cell at 50 M cells).
* The default `auto` policy keeps large runs near the no-pairwise cost by switching
  to a seeded sample of 20,000 pairs above 1,000 sequences.
* Exact all-vs-all identity for 5,000 sequences (12.5 M pairs) costs ~4 s extra
  via blocked matrix products; memory is dominated by the N × N matrices.

Reproduce:

```bash
python benchmarks/run_benchmarks.py --sizes 100,1000,5000 --lengths 1000,10000 \
    --modes stats_none,stats_auto,stats_exact,profile_basic --out benchmarks/results.csv
```

Still to do for the Application Note: comparison of shared metrics and runtime
with AMAS, AliStat and CIAlign on public alignments.
