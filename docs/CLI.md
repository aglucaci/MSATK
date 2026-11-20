# MSATK CLI

All commands share one analysis engine, so the CLI, the Python API and the reports
always agree. Run `msatk <command> --help` for every option.

| Command | Purpose |
|---|---|
| `msatk profile` | Full profile: statistics, QC, figures, offline HTML report |
| `msatk stats` | One statistics table (summary, site, sequence, pairwise, composition, gaps) to stdout/file |
| `msatk qc` | QC findings and flag tables |
| `msatk plot` | Residue-level alignment figure (PNG/SVG/PDF/TIFF) |
| `msatk logo` | Sequence logo for a region |
| `msatk encode` | ML tensors (NPZ) + manifest |
| `msatk image` | RGB or lossless label image + manifest |
| `msatk embed` | 2-D sequence embedding (PCA, UMAP, t-SNE) |
| `msatk codon` / `msatk protein` | Molecule-specific profiles |
| `msatk report` | Rebuild the report from an output directory |
| `msatk metrics` | Print the metric dictionary |
| `msatk demo` | Profile a bundled example |

Common input options: `--type {auto,dna,rna,protein,codon,translated_cds,mixed}`,
`--format auto`, `--genetic-code N` (NCBI table, default 1) and `--partitions FILE`
(RAxML/IQ-TREE partition file or NEXUS charsets). Inputs may be gzip-compressed. See
[Alignment types](alignment_types.md) and [Compatibility](compatibility.md) for every
supported format and data model.

## inspect and formats

```bash
msatk inspect genomes.maf.gz          # format, data model, contents, valid analyses
msatk inspect reads.bam --json
msatk formats                         # compatibility matrix (Markdown)
msatk formats --format json           # machine-readable capability registry
```

`inspect` reports the detected format and support level, the data model, what the file
contains, the specialized analyses that apply (and why others do not), and what is not
computed by design.

## profile

```bash
msatk profile alignment.fasta --out results/
msatk profile alignment.fasta --preset basic            # fastest: no pairwise comparison
msatk profile alignment.fasta --preset comprehensive    # exact pairwise, PCA, SVG, ML tensors
msatk profile alignment.fasta --plots png,svg,pdf --self-contained
msatk profile alignment.fasta --pairwise sampled --sample-pairs 50000 --seed 1
msatk profile alignment.fasta --embeddings umap --encodings onehot,features
```

| Option | Meaning |
|---|---|
| `--preset basic\|standard\|comprehensive` | Computation profile (see below) |
| `--pairwise auto\|exact\|sampled\|none` | Pairwise identity policy. `auto` = exact up to `--max-exact-sequences` (1000), deterministic sample of `--sample-pairs` (20,000) above |
| `--plots png,svg,pdf,tiff\|none` | Figure formats. Unsupported formats are rejected |
| `--no-alignment-figure` | Skip the residue-level overview and logo |
| `--embeddings none\|pca\|umap\|tsne` | Optional embedding (off by default). Missing dependencies raise an error rather than silently substituting PCA |
| `--encodings categorical,onehot,...` | Also export ML tensors to `ml/` |
| `--self-contained` | Embed figures so `report.html` is a single portable file |
| `--strict` | Fail on duplicate IDs, unequal lengths or invalid characters |
| `--max-gap-seq`, `--max-gap-site` | QC thresholds |
| `--genetic-code N` | NCBI translation table for stop codons, frames, codon usage and RSCU |
| `--partitions FILE` | Per-partition statistics (`partitions.csv`) |
| `--reference FASTA` | Reference for CRAM input |

`profile` dispatches on the data model. Sequence alignments get the full profile below.
Block (MAF/XMFA/AXT), read (SAM/BAM/CRAM), mapping (PAF/PSL/BLAST/GAF), graph (GFA) and
profile-HMM inputs get `summary.json`, `tables/*.csv`, `report.md`, `report.html` and
`run_manifest.json` with statistics specific to that model.

Presets:

| Preset | Pairwise | Figures | Extras |
|---|---|---|---|
| `basic` | none | PNG, no alignment overview | — |
| `standard` (default behaviour) | auto | PNG + alignment overview + logo | — |
| `comprehensive` | exact | PNG + SVG | PCA embedding, categorical + one-hot tensors |

## stats

```bash
msatk stats alignment.fasta --level site > sites.csv
msatk stats alignment.fasta --level sequence --output-format tsv --out sequences.tsv
msatk stats alignment.fasta --level summary --output-format json
msatk stats alignment.fasta --level pairwise --out pairs.csv
```

## qc

```bash
msatk qc alignment.fasta --out qc/ --max-gap-seq 0.3 --max-gap-site 0.5
```

Prints severity-ranked findings and writes `qc_findings.{csv,json}`, `qc_flags.csv` and `qc_warnings.txt`.

## plot and logo

```bash
msatk plot alignment.fasta --start 1 --end 300 --out region.svg
msatk plot proteins.faa --order consensus_identity --highlight-differences --tracks conservation,entropy,gaps --out fig.pdf
msatk plot alignment.fasta --theme my_colours.yaml --no-letters --max-rows 100 --out overview.png
msatk logo alignment.fasta --start 40 --end 80 --out logo.pdf
```

Themes: `nucleotide`, `nucleotide_cvd` (colour-vision-deficiency safe), `protein`,
`grayscale`, or a JSON/YAML file mapping symbols to hex colours.

## encode and image

```bash
msatk encode alignment.fasta --encoding onehot --out ml/
msatk encode alignment.fasta --encoding onehot --ambiguity fractional --window 256 --stride 128 --out ml/
msatk encode proteins.faa --encoding physicochemical --label family=PF00001 --out ml/
msatk image alignment.fasta --kind rgb --scale 4 --out alignment.png
msatk image alignment.fasta --kind label --out alignment_labels.png
msatk image alignment.fasta --size 224x224 --out alignment_224.png   # lossy, recorded in manifest
```

See [ML encodings](ml_encodings.md).

## report, demo, metrics

```bash
msatk report results/ --self-contained
msatk demo
msatk metrics --format json > metrics.json
```
