# 7. Beyond multiple sequence alignments

**You will:** use the same entry point to analyse whole-genome block alignments, mapped
reads, mapping and hit tables, pangenome graphs and profile HMMs. Each is loaded into
its own data model and gets statistics that make sense for it.

**Data:** small synthetic examples of each type, bundled as `toy_*`. Each is small
enough to check every number by hand.

## One entry point

```python
from msatk import load_alignment, datasets

for name in ["toy_genomes", "toy_reads", "toy_mappings", "toy_pangenome", "toy_profile",
             "ras_pf00071"]:
    obj = load_alignment(datasets.path(name))
    print(f"{name:15} -> {type(obj).__name__:18} kind={obj.kind}")
```

`msatk inspect FILE` gives the same information on the command line, together with the
analyses that are valid and those that are deliberately left out. `msatk formats` lists
every supported format with its support level (see [Compatibility](../compatibility.md)).

## Whole-genome blocks (MAF, XMFA, AXT)

Genome aligners produce *blocks*: local alignments with coordinates in each genome.
MSATK analyses the blocks and the genomes, and never joins the blocks into one matrix.

```python
maf = load_alignment(datasets.path("toy_genomes"))
s = maf.summary()
for key in ["blocks", "genomes", "reference_genome", "n50_block_columns",
            "blocks_with_all_genomes", "reference_aligned_fraction", "minus_strand_rows",
            "collinearity_breakpoints"]:
    print(f"{key:28} {s[key]}")
for g in maf.genome_rows():
    print({k: g[k] for k in ("genome", "blocks_present", "covered_bases", "aligned_fraction",
                             "minus_strand_rows", "collinearity_breakpoints")})
```

The third block is on the minus strand in `mm10`, which shows up as a strand change and
one breakpoint in collinearity with the reference. MAF minus-strand coordinates are
converted to forward-strand intervals.

## Mapped reads (SAM, BAM, CRAM)

Mapped reads get mapping QC, not alignment-column statistics:

```python
reads = load_alignment(datasets.path("toy_reads"))
s = reads.summary()
for key in ["records", "primary_reads", "mapping_rate", "mean_mapq", "mapq_0_fraction",
            "indel_rate", "mismatch_rate", "clipped_read_fraction", "spliced_reads",
            "median_insert_size", "breadth_1x"]:
    print(f"{key:24} {s[key]}")
for row in reads.reference_rows():
    print(row)
```

Depth counts aligned bases (`M`, `=` and `X` in the CIGAR string) of primary,
non-duplicate reads, as `samtools depth` does by default. The mismatch rate comes from
the `NM` tag minus indels. Reading BAM and CRAM needs `pysam`; for CRAM, pass the
reference with `--reference`.

## Mapping and hit tables (PAF, PSL, BLAST tabular, GAF)

```python
paf = load_alignment(datasets.path("toy_mappings"))
s = paf.summary()
for key in ["records", "queries", "multi_hit_queries", "mean_identity",
            "length_weighted_identity", "mean_mapq", "mean_query_coverage"]:
    print(f"{key:26} {s[key]}")
print(paf.tables()["queries"])
```

Identity follows each format's own definition: matches / block length for PAF and GAF,
`pident` for BLAST, and (matches + repMatches) / aligned bases for PSL.

## Pangenome graphs (GFA 1, rGFA)

A variation graph is analysed as a graph. Segments are classed as core, shell or private
according to how many haplotype paths pass through them.

```python
gfa = load_alignment(datasets.path("toy_pangenome"))
s = gfa.summary()
for key in ["segments", "links", "connected_components", "tip_segments",
            "branching_segment_ends", "paths", "core_segments", "shell_segments",
            "private_segments", "core_bases_fraction", "reference_rank0_segments"]:
    print(f"{key:26} {s[key]}")
for row in gfa.path_rows():
    print(row)
```

Segment `s2` (on `hap1` only) and segment `s3` (on `hap2` and the `sampleC` walk) form a
bubble between `s1` and `s4`.

## Profile HMMs (HMMER3, HHM)

```python
hmm = load_alignment(datasets.path("toy_profile"))
print(hmm.summary())
for row in hmm.state_rows():
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()})
```

## Same outputs from the command line

`msatk profile` chooses the analysis from the data model. For these inputs it writes
`summary.json`, `tables/*.csv`, `report.md`, a small offline `report.html` and
`run_manifest.json`.

```bash
msatk datasets toy_reads toy_genomes --out data
msatk profile data/toy_reads.sam --out reads_qc
msatk profile data/toy_genomes.maf --out genome_blocks
```

## What is not supported yet

PDB/mmCIF structures, GBZ and ODGI graphs, HAL and chain files, GFA 2 and BLAST XML are
*planned*. `msatk inspect` recognises them and explains how to convert them. Unaligned
FASTA/FASTQ, VCF, BED and GFF/GTF are not alignments, so MSATK does not treat them as
such.
