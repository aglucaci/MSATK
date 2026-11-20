# Benchmarks

`run_benchmarks.py` measures wall-clock time and peak resident memory of MSATK on
synthetic DNA alignments across sequence counts, alignment lengths and modes. Each
configuration runs in a separate subprocess. Results are written as CSV and JSON
(with interpreter/platform metadata). See `docs/benchmarks.md` for the latest table.
