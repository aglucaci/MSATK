# Installation

## With Conda

The recommended bioinformatics-native install path is Bioconda:

```bash
mamba install -c bioconda -c conda-forge msatk
```

or:

```bash
conda install -c bioconda -c conda-forge msatk
```

Before the official Bioconda package exists, use Conda for dependencies and pip for MSATK:

```bash
conda env create -f environment.yml
conda activate msatk
```

After Bioconda packaging is live:

```bash
conda env create -f environment-bioconda.yml
conda activate msatk
```

## With pip

```bash
pip install msatk
```

Recommended full install:

```bash
pip install "msatk[all]"
```

BAM/CRAM support requires `pysam`:

```bash
pip install "msatk[ngs]"
```

or:

```bash
mamba install -c conda-forge pysam
```

Developer install:

```bash
git clone https://github.com/yourname/msatk
cd msatk
pip install -e ".[dev,all]"
```
