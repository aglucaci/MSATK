# Quickstart

```bash
msatk profile alignment.fasta
```

MSATK automatically detects the input format and molecule type, writes tables, creates figures when plotting dependencies are available, and builds a report.

```text
alignment_msatk/
├── report.html
├── summary.json
├── parameters.yaml
├── msatk.log
├── qc_warnings.txt
├── tables/
└── figures/
```

Try the bundled example:

```bash
msatk demo
```
