# 3. Quality control of a viral coding alignment

**You will:** check the reading frame of a codon alignment, look at codon-position
variability and codon usage, use QC findings to find sequences that don't belong, and see
how removing them changes population-genetic estimates. You will also switch the genetic
code.

**Data:** 400 influenza A(H3N2) haemagglutinin coding sequences from GenBank
(1,701 nt), with a table of host labels parsed from the record titles.

## Load as a codon alignment

```python
import csv, gzip
from collections import Counter
from msatk import MSATK, datasets

ha = MSATK(datasets.path("h3n2_ha"), sequence_type="codon")  # pairwise="auto": exact for 400
s = ha.summary()
print(s["number_of_sequences"], "sequences,", s["alignment_length"] // 3, "codons")
for key in ["length_divisible_by_three", "appears_codon_aware", "stop_codon_count",
            "internal_stop_codons_exist", "frameshift_warning_count", "genetic_code"]:
    print(f"{key:28} {s[key]}")
```

Every sequence ends with a stop codon (400 terminal stops), and no internal stops or
frame problems are found, so the alignment is in frame.

## Where is the variation?

```python
from msatk.codon.analysis import codon_position_entropy, rscu

for row in codon_position_entropy(ha.stats().sites):
    print(f"codon position {row['codon_position']}: mean entropy {row['mean_entropy']:.3f} bits")
```

Third positions vary about 3.5 times more than first and second positions, as expected
when most substitutions are synonymous.

```python
leu = [r for r in rscu(ha.alignment) if r["amino_acid"] == "L"]
for r in sorted(leu, key=lambda r: -r["rscu"]):
    print(f"{r['codon']}  count {r['count']:>5}  RSCU {r['rscu']:.2f}")
```

RSCU (relative synonymous codon usage) compares each codon's count with the count
expected if all synonymous codons were used equally. 1.0 means no bias.

## What QC finds

```python
for f in ha.qc_findings():
    print(f"[{f.severity}] {f.code}: {f.message}")
```

Which sequences are divergent? The record titles give each sequence's host:

```python
with gzip.open(datasets.path("h3n2_ha_hosts"), "rt") as handle:
    host = {r["accession"]: r["host_label"] for r in csv.DictReader(handle, delimiter="\t")}

divergent = set(next(f for f in ha.qc_findings() if f.code == "DIVERGENT_SEQUENCES").items)
print("all sequences:       ", dict(Counter(host.values()).most_common()))
print("flagged as divergent:", dict(Counter(host[a] for a in divergent).most_common()))
```

QC flags all the canine, avian and environmental sequences and most of the swine ones.
Those come from other host lineages that were caught by a broad database query. The 11
flagged human sequences sit at the low end of the human identity distribution and are
worth a look, but they are not from a different host. The
swine sequences that are *not* flagged have human-like haemagglutinin genes:

```python
rows = {r["sequence_id"]: r for r in ha.stats().sequences}
for acc in sorted(a for a in host if host[a] == "swine" and a not in divergent):
    r = rows[acc]
    print(acc, f"identity to consensus {r['consensus_identity']:.3f};",
          f"nearest neighbour {r['nearest_neighbor']} ({host[r['nearest_neighbor']]},",
          f"{r['nearest_neighbor_identity']:.4f})")
```

## Does it matter?

Compare the population-genetic summaries before and after removing the flagged
sequences:

```python
from msatk.encoding.transforms import subset_alignment
from msatk.stats import compute_statistics

keep = [a for a in ha.alignment.ids if a not in divergent]
filtered, coords = subset_alignment(ha.alignment, ids=keep)
after = compute_statistics(filtered, "codon", pairwise="none").summary

print(f"{'':26}{'all 400':>10}{'after QC':>10}")
for key in ["variable_sites", "nucleotide_diversity_pi", "theta_watterson", "tajimas_d",
            "ts_tv_ratio"]:
    print(f"{key:26}{s[key]:>10.4g}{after[key]:>10.4g}")
print("rows kept:", len(coords.original_rows))
```

Removing 27 sequences halves nucleotide diversity and Watterson's θ. Any later
analysis that depends on diversity, such as selection tests or dating, would otherwise
be inflated by mixing hosts. `subset_alignment` also returns a coordinate map, so the
filtered alignment can be traced back to the original rows and sites.

## Other genetic codes

All NCBI translation tables are available, and they affect stop codons, frame QC, codon
usage and RSCU. In the vertebrate mitochondrial code (table 2), TGA codes for tryptophan:

```python
from pathlib import Path
Path("mito.fasta").write_text(">a\nATGTGAAAATAA\n>b\nATGTGGAAATAA\n")

for code in (1, 2):
    m = MSATK("mito.fasta", sequence_type="codon", genetic_code=code)
    stops = [f.message for f in m.qc_findings() if f.code == "INTERNAL_STOP_CODON"]
    print(f"genetic code {code}: {stops or 'no internal stop codons'}")
```

On the command line, use `msatk profile cds.fasta --type codon --genetic-code 2`.
`msatk.core.genetic_code.available_codes()` lists every table.
