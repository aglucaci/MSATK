# Python API

```python
from msatk import AlignmentProfiler, profile_alignment

profiler = AlignmentProfiler("alignment.fasta")
summary = profiler.summary()
per_sequence = profiler.per_sequence_stats()
per_site = profiler.per_site_stats()

profiler.plot_gap_profile()
profiler.plot_entropy()
profiler.write_report("report.html")

results = profile_alignment("alignment.fasta", outdir="alignment_msatk")
```

When pandas is installed, table-like methods return pandas DataFrames. In minimal environments, they return lists of dictionaries.
