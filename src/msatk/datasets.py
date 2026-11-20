"""Small bundled datasets for tutorials and testing.

>>> from msatk import datasets
>>> datasets.list_datasets()
>>> path = datasets.path("ras_pf00071")

Real datasets are public alignments redistributed unchanged (gzip-compressed);
``toy_*`` datasets are small synthetic files that illustrate each alignment type.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib import resources
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Dataset:
    name: str
    filename: str
    kind: str
    description: str
    source: str
    license: str


DATASETS: dict[str, Dataset] = {
    d.name: d
    for d in (
        Dataset(
            "ras_pf00071",
            "ras_pf00071.sto.gz",
            "protein MSA",
            "Ras family seed alignment, 60 sequences x 229 columns (Stockholm).",
            "Pfam PF00071 seed via InterPro (Mistry et al. 2021)",
            "CC0",
        ),
        Dataset(
            "h3n2_ha",
            "h3n2_ha.fasta.gz",
            "codon MSA",
            "Influenza A(H3N2) haemagglutinin CDS, 400 x 1,701 nt (GenBank records "
            "published in 2017; includes swine, canine, avian and environmental hosts).",
            "NCBI GenBank (nuccore), see h3n2_ha_hosts for accessions and titles",
            "NCBI: no restrictions on use",
        ),
        Dataset(
            "h3n2_ha_hosts",
            "h3n2_ha_hosts.tsv.gz",
            "metadata",
            "Accession, host label (parsed from the strain name) and GenBank title for h3n2_ha.",
            "NCBI GenBank record titles",
            "NCBI: no restrictions on use",
        ),
        Dataset(
            "trna_rf00005",
            "trna_rf00005.sto.gz",
            "RNA MSA with structure",
            "tRNA seed alignment, 954 sequences x 118 columns with consensus "
            "secondary structure (#=GC SS_cons).",
            "Rfam RF00005 seed",
            "CC0",
        ),
        Dataset("toy_genomes", "toy_genomes.maf", "block alignment",
                "Three genomes, three MAF blocks, one inverted block.", "synthetic", "CC0"),
        Dataset("toy_reads", "toy_reads.sam", "read alignment",
                "Eight SAM records: pairs, clipping, indels, splicing, duplicate, secondary, "
                "unmapped.", "synthetic", "CC0"),
        Dataset("toy_mappings", "toy_mappings.paf", "mapping",
                "Three PAF records, one multi-mapping read.", "synthetic", "CC0"),
        Dataset("toy_pangenome", "toy_pangenome.gfa", "graph",
                "rGFA graph with a bubble, two paths and one walk.", "synthetic", "CC0"),
        Dataset("toy_profile", "toy_profile.hmm", "profile HMM",
                "Three-state HMMER3 DNA profile.", "synthetic", "CC0"),
        Dataset("toy_pairwise", "toy_pairwise.needle", "pairwise alignment",
                "EMBOSS needle global protein alignment.", "synthetic", "CC0"),
        Dataset("toy_partitioned", "toy_partitioned.nexus", "partitioned MSA",
                "NEXUS DNA matrix with gene and codon-position charsets.", "synthetic", "CC0"),
    )
}  # fmt: skip


def list_datasets() -> list[dict[str, Any]]:
    """Describe every bundled dataset."""

    return [asdict(d) for d in DATASETS.values()]


def path(name: str) -> Path:
    """Filesystem path of a bundled dataset (``.gz`` files can be read directly)."""

    try:
        filename = DATASETS[name].filename
    except KeyError as exc:
        raise KeyError(f"Unknown dataset '{name}'. Available: {', '.join(DATASETS)}") from exc
    return Path(str(resources.files("msatk.data.datasets").joinpath(filename)))
