"""Machine-readable capability registry for alignment formats and data models.

Every format MSATK knows about is declared here with its parser, data model,
alphabets, preserved annotations, compatible analyses and a support level:

``stable``   implemented, covered by tests on real or fixture data, P0 scope
``beta``     implemented and tested on small fixtures; interfaces may change
``experimental`` implemented but not yet covered by a test fixture
``planned``  recognised and reported by ``msatk inspect``, but not yet readable

``docs/compatibility.md`` is generated from this registry (``msatk formats``) and a
test keeps the two in sync, so documentation cannot claim unvalidated support.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

SEQUENCE_MODULES = (
    "statistics",
    "qc",
    "figures",
    "html_report",
    "ml_encodings",
)

#: Data models (alignment kinds) and the analyses valid for each.
KINDS: dict[str, dict[str, Any]] = {
    "sequence": {
        "model": "SequenceAlignment",
        "description": "Rows of one rectangular alignment (pairwise or multiple).",
        "analyses": [
            "alignment/site/sequence/pair statistics",
            "QC findings",
            "population genetics (nucleotides)",
            "codon usage, RSCU, frames, alternative genetic codes (codon)",
            "physicochemical and BLOSUM62 substitution scores (protein)",
            "pairwise-alignment statistics (exactly two sequences)",
            "RNA base-pair conservation and covariation (SS_cons)",
            "profile match/insert states (A2M/A3M, #=GC RF)",
            "per-partition statistics (charsets / partition file)",
            "alignment figures, logos, HTML report, ML encodings",
        ],
        "not_valid": [],
    },
    "block": {
        "model": "BlockAlignment",
        "description": "Local alignment blocks with genomic coordinates (whole-genome).",
        "analyses": [
            "block size distribution (N50)",
            "per-genome aligned fraction, coverage and duplication",
            "strand and collinearity breakpoints relative to the reference",
            "per-block gap, variability and consensus identity",
        ],
        "not_valid": [
            "concatenated MSA statistics (blocks are never merged)",
            "all-vs-all identity across genomes",
        ],
    },
    "read": {
        "model": "ReadAlignment",
        "description": "Sequencing reads mapped to reference coordinates (SAM/BAM/CRAM).",
        "analyses": [
            "mapping rate and flag summary",
            "MAPQ distribution",
            "CIGAR: clipping, indels, splicing, mismatch rate",
            "base quality, read length, insert size",
            "per-reference depth and breadth of coverage",
        ],
        "not_valid": ["all-vs-all read identity", "column entropy over a read pileup"],
    },
    "mapping": {
        "model": "MappingAlignment",
        "description": "Pairwise local hits/mappings (PAF, PSL, BLAST tabular, GAF).",
        "analyses": [
            "hit counts and multi-hit queries",
            "identity (format-specific definition) and alignment length",
            "MAPQ, e-value",
            "query and target coverage",
        ],
        "not_valid": ["column-wise MSA statistics"],
    },
    "graph": {
        "model": "GraphAlignment",
        "description": "Sequence/variation graphs and pangenome paths (GFA 1, rGFA).",
        "analyses": [
            "segment, link and component statistics",
            "tips and branching",
            "path/walk lengths",
            "core / shell / private segments across paths",
        ],
        "not_valid": ["flattening into a rectangular alignment"],
    },
    "profile_model": {
        "model": "ProfileModel",
        "description": "Profile HMMs (HMMER3, HH-suite HHM).",
        "analyses": ["model metadata", "match-state emission entropy (HMMER3)"],
        "not_valid": ["sequence-level statistics (no sequences in the file)"],
    },
    "structure": {
        "model": "StructureAlignment (planned)",
        "description": "3D coordinates and structural superpositions (PDB, mmCIF).",
        "analyses": [],
        "not_valid": ["all analyses until structural alignment support is implemented"],
    },
}


@dataclass(frozen=True)
class FormatCapability:
    name: str
    label: str
    extensions: tuple[str, ...]
    kind: str
    level: str
    reader: str
    alphabets: tuple[str, ...] = ("dna", "rna", "protein")
    preserves: tuple[str, ...] = ()
    writer: bool = False
    gzip: bool = True
    dependency: str | None = None
    tested_with: tuple[str, ...] = ()
    notes: str = ""
    family: str = ""
    extra: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["model"] = KINDS[self.kind]["model"]
        return out


_SEQ = ("dna", "rna", "protein")
F = FormatCapability
FORMATS: tuple[FormatCapability, ...] = (
    # ---------------------------------------------------------------- standard MSA
    F(
        "fasta",
        "FASTA (aligned)",
        ("fa", "fasta", "fas", "faa", "fna", "afa", "mfa"),
        "sequence",
        "stable",
        "native",
        _SEQ,
        ("ids", "descriptions"),
        family="Standard MSA",
        tested_with=("tests/data/formats/tiny.aligned.fasta",),
    ),
    F(
        "phylip",
        "PHYLIP (strict and relaxed, sequential)",
        ("phy", "phylip"),
        "sequence",
        "stable",
        "native",
        _SEQ,
        ("ids",),
        family="Standard MSA",
        tested_with=("tests/data/formats/tiny.phy", "tests/data/formats/tiny_relaxed.phylip"),
    ),
    F(
        "clustal",
        "CLUSTAL",
        ("aln", "clustal", "clw"),
        "sequence",
        "stable",
        "native",
        _SEQ,
        ("ids",),
        family="Standard MSA",
        tested_with=("tests/data/formats/tiny.clustal",),
    ),
    F(
        "nexus",
        "NEXUS",
        ("nex", "nexus", "nxs"),
        "sequence",
        "stable",
        "native",
        _SEQ,
        ("ids", "charset partitions"),
        family="Standard MSA / multi-locus",
        tested_with=("tests/data/formats/tiny.nexus", "tests/data/formats/partitioned.nexus"),
    ),
    F(
        "stockholm",
        "Stockholm",
        ("sto", "stk", "stockholm"),
        "sequence",
        "stable",
        "native",
        _SEQ,
        ("ids", "#=GC SS_cons", "#=GC RF", "other #=GC tracks"),
        family="Standard MSA / RNA structure",
        tested_with=("tests/data/formats/tiny.stockholm", "tests/data/formats/rna_structure.sto"),
    ),
    # ---------------------------------------------------------------- specialized MSA
    F(
        "a3m",
        "A3M (HH-suite)",
        ("a3m",),
        "sequence",
        "beta",
        "native",
        ("protein", "dna", "rna"),
        ("match states", "insert lengths"),
        family="Specialized MSA / profile",
        tested_with=("tests/data/formats/tiny.a3m",),
        notes="Lowercase = insertions (removed from the matrix, counted); insert gaps omitted.",
    ),
    F(
        "a2m",
        "A2M",
        ("a2m",),
        "sequence",
        "beta",
        "native",
        ("protein", "dna", "rna"),
        ("match states", "insert lengths"),
        family="Specialized MSA / profile",
        tested_with=("tests/data/formats/tiny.a2m",),
        notes="Lowercase = insertions, '.' = insert gaps.",
    ),
    F(
        "msf",
        "MSF / GCG",
        ("msf",),
        "sequence",
        "beta",
        "native",
        _SEQ,
        ("ids",),
        family="Specialized MSA",
        tested_with=("tests/data/formats/tiny.msf",),
        notes="'.' and '~' are read as gaps.",
    ),
    F(
        "mega",
        "MEGA",
        ("meg", "mega"),
        "sequence",
        "beta",
        "native",
        _SEQ,
        ("ids",),
        family="Specialized MSA",
        tested_with=("tests/data/formats/tiny.meg",),
        notes="'.' (identical to first sequence) is expanded.",
    ),
    F(
        "selex",
        "SELEX",
        ("slx", "selex"),
        "sequence",
        "beta",
        "native",
        _SEQ,
        ("ids",),
        family="Specialized MSA",
        tested_with=("tests/data/formats/tiny.slx",),
    ),
    F(
        "mase",
        "MASE",
        ("mase",),
        "sequence",
        "beta",
        "native",
        _SEQ,
        ("ids",),
        family="Specialized MSA",
        tested_with=("tests/data/formats/tiny.mase",),
    ),
    # ---------------------------------------------------------------- pairwise
    F(
        "emboss",
        "EMBOSS pairwise (needle, water, stretcher, matcher)",
        ("needle", "water", "emboss"),
        "sequence",
        "beta",
        "native",
        _SEQ,
        ("method (global/local)", "reported score", "scoring parameters"),
        family="Pairwise alignment",
        tested_with=("tests/data/formats/tiny.needle",),
    ),
    F(
        "paf",
        "PAF",
        ("paf",),
        "mapping",
        "beta",
        "native",
        ("dna", "rna"),
        ("coordinates", "strand", "MAPQ"),
        family="Pairwise alignment / mapping",
        tested_with=("tests/data/formats/tiny.paf",),
    ),
    F(
        "psl",
        "PSL (BLAT)",
        ("psl",),
        "mapping",
        "beta",
        "native",
        _SEQ,
        ("coordinates", "strand", "block counts"),
        family="Pairwise alignment / mapping",
        tested_with=("tests/data/formats/tiny.psl",),
    ),
    F(
        "blast_tab",
        "BLAST/DIAMOND tabular (-outfmt 6)",
        ("blast6", "outfmt6", "m8"),
        "mapping",
        "beta",
        "native",
        _SEQ,
        ("coordinates", "e-value", "bit score"),
        family="Pairwise alignment / hits",
        tested_with=("tests/data/formats/tiny.m8",),
    ),
    F(
        "blast_xml",
        "BLAST XML / JSON",
        ("xml", "json"),
        "mapping",
        "planned",
        "none",
        _SEQ,
        family="Pairwise alignment / hits",
        gzip=False,
        notes="Use -outfmt 6 for now.",
    ),
    # ---------------------------------------------------------------- reads
    F(
        "sam",
        "SAM",
        ("sam",),
        "read",
        "beta",
        "native",
        ("dna", "rna"),
        ("reference coordinates", "CIGAR", "flags", "MAPQ", "NM", "qualities"),
        family="Sequencing reads",
        tested_with=("tests/data/formats/tiny.sam", "tests/data/formats/reads.sam"),
    ),
    F(
        "bam",
        "BAM",
        ("bam",),
        "read",
        "beta",
        "pysam",
        ("dna", "rna"),
        ("reference coordinates", "CIGAR", "flags", "MAPQ", "NM", "qualities"),
        family="Sequencing reads",
        gzip=False,
        dependency="pysam",
        tested_with=("tests/data/formats/tiny.bam",),
    ),
    F(
        "cram",
        "CRAM",
        ("cram",),
        "read",
        "experimental",
        "pysam",
        ("dna", "rna"),
        ("reference coordinates", "CIGAR", "flags", "MAPQ", "NM", "qualities"),
        family="Sequencing reads",
        gzip=False,
        dependency="pysam",
        notes="Pass the reference FASTA when the CRAM needs it.",
    ),
    # ---------------------------------------------------------------- whole genome
    F(
        "maf",
        "MAF (multiple alignment format)",
        ("maf",),
        "block",
        "beta",
        "native",
        ("dna",),
        ("block membership", "coordinates", "strand", "source sizes", "scores"),
        family="Whole-genome alignment",
        tested_with=("tests/data/formats/tiny.maf", "tests/data/formats/blocks.maf"),
    ),
    F(
        "xmfa",
        "XMFA (Mauve)",
        ("xmfa",),
        "block",
        "beta",
        "native",
        ("dna",),
        ("block membership", "coordinates", "strand"),
        family="Whole-genome alignment",
        tested_with=("tests/data/formats/tiny.xmfa",),
    ),
    F(
        "axt",
        "AXT",
        ("axt",),
        "block",
        "beta",
        "native",
        ("dna",),
        ("coordinates", "strand", "scores"),
        family="Whole-genome alignment",
        tested_with=("tests/data/formats/tiny.axt",),
    ),
    F(
        "chain",
        "UCSC chain",
        ("chain",),
        "block",
        "planned",
        "none",
        ("dna",),
        family="Whole-genome alignment",
        notes="No sequences; needs a block-size model.",
    ),
    F(
        "hal",
        "HAL (Cactus)",
        ("hal",),
        "block",
        "planned",
        "none",
        ("dna",),
        gzip=False,
        family="Whole-genome alignment",
        notes="Export to MAF with hal2maf.",
    ),
    # ---------------------------------------------------------------- profiles
    F(
        "hmmer",
        "HMMER3 profile HMM",
        ("hmm",),
        "profile_model",
        "beta",
        "native",
        ("protein", "dna"),
        ("model metadata", "match emissions", "COMPO"),
        family="Profile alignment",
        tested_with=("tests/data/formats/tiny.hmm",),
    ),
    F(
        "hhm",
        "HH-suite HHM",
        ("hhm",),
        "profile_model",
        "beta",
        "native",
        ("protein",),
        ("model metadata",),
        family="Profile alignment",
        tested_with=("tests/data/formats/tiny.hhm",),
        notes="Metadata only.",
    ),
    # ---------------------------------------------------------------- graphs
    F(
        "gfa",
        "GFA 1 / rGFA",
        ("gfa", "gfa1", "rgfa"),
        "graph",
        "beta",
        "native",
        ("dna",),
        ("segments", "links", "paths", "walks", "rGFA rank"),
        family="Pangenome graph",
        tested_with=("tests/data/formats/tiny.gfa",),
        notes="GFA 2 not yet supported.",
    ),
    F(
        "gaf",
        "GAF (graph alignment)",
        ("gaf",),
        "mapping",
        "beta",
        "native",
        ("dna",),
        ("graph paths", "coordinates", "MAPQ"),
        family="Pangenome graph",
        tested_with=("tests/data/formats/tiny.gaf",),
    ),
    F(
        "gbz",
        "GBZ",
        ("gbz",),
        "graph",
        "planned",
        "none",
        ("dna",),
        gzip=False,
        family="Pangenome graph",
        notes="Convert with `vg view` / `gbz2gfa`.",
    ),
    F(
        "odgi",
        "ODGI",
        ("og",),
        "graph",
        "planned",
        "none",
        ("dna",),
        gzip=False,
        family="Pangenome graph",
        notes="Convert with `odgi view -g`.",
    ),
    # ---------------------------------------------------------------- structure
    F(
        "pdb",
        "PDB",
        ("pdb", "ent"),
        "structure",
        "planned",
        "none",
        ("protein",),
        family="Protein structure",
        notes="Coordinates are not a residue alignment.",
    ),
    F(
        "mmcif",
        "mmCIF",
        ("cif", "mmcif"),
        "structure",
        "planned",
        "none",
        ("protein",),
        family="Protein structure",
        notes="Coordinates are not a residue alignment.",
    ),
)

_BY_NAME = {f.name: f for f in FORMATS}


def get_format(name: str) -> FormatCapability:
    try:
        return _BY_NAME[name]
    except KeyError as exc:
        raise KeyError(f"Unknown format '{name}'. Known: {', '.join(sorted(_BY_NAME))}") from exc


def capability_registry() -> dict[str, Any]:
    return {
        "kinds": KINDS,
        "formats": [f.as_dict() for f in FORMATS],
        "levels": {
            "stable": "Implemented and tested; part of the publication-ready scope.",
            "beta": "Implemented and tested on fixtures; interface may change.",
            "experimental": "Implemented but not yet covered by a test fixture.",
            "planned": "Recognised by `msatk inspect`, not yet readable.",
        },
    }


def compatibility_markdown() -> str:
    """Render the registry as Markdown (``docs/compatibility.md``)."""

    lines = [
        "# Format and alignment-type compatibility",
        "",
        "<!-- Generated by `msatk formats > docs/compatibility.md`; do not edit by hand. -->",
        "",
        "MSATK distinguishes the **alignment type** (data model), the **method** that produced",
        "it (global, local, ...; recorded as metadata when the format states it) and the",
        "**file format**. Each format below is declared in a machine-readable registry",
        "(`msatk.alignments.registry`), and this page is generated from it.",
        "",
        "Support levels: **stable** (tested, publication scope), **beta** (implemented and",
        "tested on fixtures), **experimental** (implemented, no test fixture yet), **planned**",
        "(recognised by `msatk inspect` but not readable yet).",
        "Gzip-compressed input (`.gz`) is accepted for every text format marked *gz*.",
        "",
        "## Data models",
        "",
        "| Kind | Model | Description | Valid analyses | Deliberately not computed |",
        "|---|---|---|---|---|",
    ]
    for kind, info in KINDS.items():
        lines.append(
            f"| `{kind}` | {info['model']} | {info['description']} | "
            f"{'; '.join(info['analyses']) or '–'} | {'; '.join(info['not_valid']) or '–'} |"
        )
    lines += [
        "",
        "## Formats",
        "",
        "| Format | Family | Extensions | Kind | Level | Reader | gz | Preserved annotations | Notes |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for f in FORMATS:
        reader = f.reader + (
            f" ({f.dependency})" if f.dependency and f.reader != f.dependency else ""
        )
        lines.append(
            f"| {f.label} | {f.family} | {', '.join('.' + e for e in f.extensions)} | `{f.kind}` | "
            f"{f.level} | {reader} | {'gz' if f.gzip else '–'} | {', '.join(f.preserves) or '–'} | "
            f"{f.notes or '–'} |"
        )
    lines += [
        "",
        "Related inputs that are **not** alignments — unaligned FASTA/FASTQ, VCF, BED, GFF/GTF —",
        "are reads, sequences, variants or annotations and are not treated as completed",
        "alignments.",
        "",
    ]
    return "\n".join(lines)
