# Quickstart

![MSATK logo](assets/msatk_logo.png)

```bash
mamba install -c bioconda -c conda-forge msatk
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

Common alignment formats are detected automatically:

```bash
msatk profile alignment.fasta
msatk profile alignment.phy
msatk profile alignment.aln
msatk profile alignment.sto
msatk profile alignment.nex
msatk profile alignment.maf
msatk profile reads.sam
msatk profile reads.bam
msatk profile reads.cram
```

BAM/CRAM support requires `pysam`:

```bash
pip install "msatk[ngs]"
```
