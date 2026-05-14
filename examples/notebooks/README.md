# MSATK Notebooks

This directory is reserved for notebook examples. Start with:

```python
from msatk import AlignmentProfiler

profiler = AlignmentProfiler("../dna/tiny_dna.fasta")
summary = profiler.summary()
per_site = profiler.per_site_stats()
```
