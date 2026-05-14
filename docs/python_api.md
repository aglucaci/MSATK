# Python API

```python
from msatk import MSATK, profile_alignment

profiler = MSATK("alignment.fasta")
summary = profiler.summary()
per_sequence = profiler.per_sequence_stats()
per_site = profiler.per_site_stats()

profiler.plot_gap_profile()
profiler.plot_entropy()
profiler.write_report("report.html")

results = profile_alignment("alignment.fasta", outdir="alignment_msatk")
```

When pandas is installed, table-like methods return pandas DataFrames. In minimal environments, they return lists of dictionaries.

Lowercase alias:

```python
from msatk import msatk

profiler = msatk("alignment.fasta")
```
