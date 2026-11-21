# Python API

```python
from msatk import MSATK, AlignmentProfiler, profile_alignment

msa = AlignmentProfiler("alignment.fasta", pairwise="auto")   # alias of MSATK

# Statistics (DataFrames when pandas is installed, else lists of dicts)
summary = msa.summary()
sites = msa.sites()
sequences = msa.sequences()
pairs = msa.pairwise()          # long table; computes the exact matrix on demand
gaps = msa.gap_blocks()
bundle = msa.stats()            # StatisticsBundle: summary, sites, sequences, pairwise, encoded

# Quality control
findings = msa.qc()             # list of dicts, severity-ranked

# Human visualisation
msa.plot_alignment("alignment.svg", start=1, end=300, theme="nucleotide")
msa.plot_logo("logo.pdf", start=1, end=60)

# Machine learning
tensor = msa.to_tensor("onehot")              # EncodedTensor (.array, .channels, .ids)
msa.to_image("alignment.png", kind="label")
msa.export_encoding("ml/", encoding="onehot", strategy="tiles", height=128, width=256)

# Reports
msa.report("report.html")                     # interactive, single file
results = profile_alignment("alignment.fasta", outdir="alignment_msatk", preset="standard")
```

Any alignment type:

```python
from msatk import load_alignment, inspect_alignment

aln = load_alignment("reads.bam")        # ReadAlignment
aln.kind, aln.summary()["mapping_rate"]
aln.tables()["references"]               # per-reference depth and breadth

inspect_alignment("family.sto")          # format, model, annotations, applicable analyses
seq = load_alignment("family.sto")       # SequenceAlignment
seq.summary()["specialized"]["rna_structure"]

MSATK("cds.fasta", sequence_type="codon", genetic_code=2, partitions="parts.txt").specialized()
```

Lower-level modules:

* `msatk.stats` — `compute_statistics`, `per_site_stats`, `per_sequence_stats`,
  `alignment_summary`, `compute_pairwise`, `diversity_summary`, `shannon_entropy`;
* `msatk.stats.registry` — the metric dictionary;
* `msatk.alignments` — data models, `load_alignment`, `inspect_alignment` and the
  capability registry;
* `msatk.stats.specialized` — pairwise, substitution, RNA-structure, profile-state and
  partition analyses;
* `msatk.core.genetic_code` — NCBI translation tables;
* `msatk.core.alphabet` — character policies per molecule type;
* `msatk.encoding` — tensors, images, transforms, export;
* `msatk.viz` — alignment figures, logos, themes;
* `msatk.ml` — embeddings and dataset adapters.

Earlier method names (`per_sequence_stats`, `per_site_stats`, `pairwise_identity_matrix`,
`distance_matrix`, `qc_warnings`, `write_report`) remain available.
