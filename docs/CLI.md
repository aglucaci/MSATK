# MSATK CLI

The primary command is:

```bash
msatk profile alignment.fasta
```

It writes a FastQC-like directory named `<input_stem>_msatk/` unless `--out` is provided.

Supported inputs include FASTA/aligned FASTA, A3M, PHYLIP/relaxed PHYLIP, CLUSTAL, Stockholm, NEXUS, MAF, SAM, BAM, and CRAM. BAM/CRAM support requires `pysam`.

MSATK automatically reports detected file format, molecule type, alignment length, length consistency, stop-codon presence, frameshift warnings, and codon-aware status in `summary.json`.

## Profile

```bash
msatk profile alignment.fasta --out msatk_results/
```

## QC

```bash
msatk qc alignment.fasta --out qc_results/ --max-gap-seq 0.3 --max-gap-site 0.5
```

## Codon

```bash
msatk codon cds_alignment.fasta --out codon_results/
```

## Protein

```bash
msatk protein proteins.faa --out protein_results/
```

## Embeddings

```bash
msatk embed alignment.fasta --method pca --representation onehot --out embeddings/
```

## Report

```bash
msatk report msatk_results/ --format html
```

## Demo

```bash
msatk demo
```
