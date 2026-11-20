# MSATK examples

Runnable, worked examples are in the tutorials:
[`docs/tutorials/`](../../docs/tutorials/index.md). They use datasets bundled with
the package:

```python
from msatk import MSATK, datasets

msa = MSATK(datasets.path("ras_pf00071"))
summary = msa.summary()
sites = msa.sites()
```
