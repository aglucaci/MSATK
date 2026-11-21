"""Input format detection and sequence-alignment readers.

``read_alignment`` returns a rectangular :class:`~msatk.models.Alignment` for
sequence-alignment formats. ``read_sequence_alignment`` additionally returns the
annotations a format carries (Stockholm ``#=GC`` tracks such as ``SS_cons`` and
``RF``, A2M/A3M insert states, NEXUS character sets, EMBOSS pairwise metadata).
Gzip-compressed input is accepted everywhere. Block, read, mapping, graph and
profile-model formats are handled by :func:`msatk.alignments.load_alignment`.
"""

from __future__ import annotations

import contextlib
import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from msatk.exceptions import AlignmentFormatError
from msatk.io.validate import validate_alignment
from msatk.models import Alignment, SequenceRecord

#: Extensions (after removing ``.gz``/``.bgz``) mapped to canonical format names.
EXTENSION_FORMATS: dict[str, str] = {
    "fa": "fasta",
    "fasta": "fasta",
    "fas": "fasta",
    "faa": "fasta",
    "fna": "fasta",
    "ffn": "fasta",
    "afa": "fasta",
    "mfa": "fasta",
    "a2m": "a2m",
    "a3m": "a3m",
    "phy": "phylip",
    "phylip": "phylip",
    "aln": "clustal",
    "clustal": "clustal",
    "clw": "clustal",
    "sto": "stockholm",
    "stk": "stockholm",
    "stockholm": "stockholm",
    "nex": "nexus",
    "nexus": "nexus",
    "nxs": "nexus",
    "msf": "msf",
    "meg": "mega",
    "mega": "mega",
    "slx": "selex",
    "selex": "selex",
    "mase": "mase",
    "needle": "emboss",
    "water": "emboss",
    "emboss": "emboss",
    "maf": "maf",
    "xmfa": "xmfa",
    "axt": "axt",
    "sam": "sam",
    "bam": "bam",
    "cram": "cram",
    "paf": "paf",
    "psl": "psl",
    "blast6": "blast_tab",
    "outfmt6": "blast_tab",
    "m8": "blast_tab",
    "gfa": "gfa",
    "gfa1": "gfa",
    "rgfa": "gfa",
    "gaf": "gaf",
    "hmm": "hmmer",
    "hhm": "hhm",
    "pdb": "pdb",
    "ent": "pdb",
    "cif": "mmcif",
    "mmcif": "mmcif",
    "gbz": "gbz",
    "og": "odgi",
    "hal": "hal",
    "chain": "chain",
}

SEQUENCE_FORMATS = {
    "fasta",
    "a2m",
    "a3m",
    "phylip",
    "clustal",
    "stockholm",
    "nexus",
    "msf",
    "mega",
    "selex",
    "mase",
    "emboss",
    "maf",
    "sam",
    "bam",
    "cram",
}


@dataclass
class Partition:
    """A named set of 1-based alignment columns (e.g. a gene or codon position)."""

    name: str
    sites: list[int]
    datatype: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {"name": self.name, "datatype": self.datatype, "n_sites": len(self.sites)}


@dataclass
class SequenceAnnotations:
    """Format-level annotations that must survive parsing."""

    gc: dict[str, str] = field(default_factory=dict)  # Stockholm #=GC tracks (SS_cons, RF, ...)
    insertions: list[list[int]] | None = None  # A2M/A3M: inserted residues after each match col
    partitions: list[Partition] = field(default_factory=list)
    pairwise: dict[str, Any] = field(default_factory=dict)  # EMBOSS metadata
    method: str = "unknown"  # global / local / unknown

    @property
    def ss_cons(self) -> str | None:
        return self.gc.get("SS_cons")

    @property
    def rf(self) -> str | None:
        return self.gc.get("RF")

    def present(self) -> list[str]:
        keys = [f"GC:{k}" for k in self.gc]
        if self.insertions is not None:
            keys.append("insert_states")
        if self.partitions:
            keys.append("partitions")
        if self.pairwise:
            keys.append("pairwise_metadata")
        return keys


def canonical_suffix(path: Path) -> str:
    """Return the format-bearing extension, ignoring ``.gz``/``.bgz`` compression."""

    suffixes = [s.lower().lstrip(".") for s in path.suffixes]
    while suffixes and suffixes[-1] in {"gz", "bgz", "gzip"}:
        suffixes.pop()
    return suffixes[-1] if suffixes else ""


def is_gzipped(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def read_text(path: str | Path) -> str:
    """Read a text file, transparently decompressing gzip/bgzip input."""

    source = Path(path)
    if is_gzipped(source):
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            return handle.read()
    return source.read_text(encoding="utf-8")


def read_alignment(
    path: str | Path, fmt: str = "auto", validation_mode: str = "permissive"
) -> Alignment:
    """Read a sequence alignment (see module docstring for formats)."""

    return read_sequence_alignment(path, fmt, validation_mode)[0]


def read_sequence_alignment(
    path: str | Path, fmt: str = "auto", validation_mode: str = "permissive"
) -> tuple[Alignment, SequenceAnnotations]:
    """Read a sequence alignment and the annotations its format carries."""

    source = Path(path)
    detected = _detect_binary_format(source) if fmt == "auto" else fmt.lower()
    detected = EXTENSION_FORMATS.get(detected, detected)
    text = ""
    ann = SequenceAnnotations()
    if detected not in {"bam", "cram"}:
        text = read_text(source)
        if not text.strip():
            raise AlignmentFormatError(f"MSATK could not read {source}: the file is empty.")
        detected = detect_format(text, source) if fmt == "auto" else detected
    if detected == "fasta":
        records = _read_fasta(text)
    elif detected in {"a2m", "a3m"}:
        records, ann.insertions = _read_a2m_a3m(text, detected)
    elif detected == "phylip":
        records = _read_phylip(text)
    elif detected == "clustal":
        records = _read_clustal(text)
    elif detected == "stockholm":
        records, ann.gc = _read_stockholm_annotated(text)
    elif detected == "nexus":
        records = _read_nexus(text)
        ann.partitions = parse_nexus_charsets(text)
    elif detected == "msf":
        records = _read_msf(text)
    elif detected == "mega":
        records = _read_mega(text)
    elif detected == "selex":
        records = _read_selex(text)
    elif detected == "mase":
        records = _read_mase(text)
    elif detected == "emboss":
        records, ann.pairwise = _read_emboss(text)
        ann.method = str(ann.pairwise.get("method", "unknown"))
    elif detected == "maf":
        records = _read_maf(text)
    elif detected == "sam":
        records = _read_sam(text)
    elif detected in {"bam", "cram"}:
        records = _read_bam_cram(source, detected)
    else:
        raise AlignmentFormatError(
            f"'{detected}' is not a sequence-alignment format. Use msatk.alignments."
            "load_alignment() or `msatk inspect` for block, read, mapping, graph and "
            "profile inputs."
        )
    alignment = Alignment(tuple(records), source=str(source), fmt=detected)
    validate_alignment(alignment, mode=validation_mode)
    return alignment, ann


def detect_format(text: str, path: Path | None = None) -> str:
    """Detect a text format from its content, falling back to the file extension."""

    stripped = text.lstrip()
    suffix = canonical_suffix(path) if path else ""
    by_ext = EXTENSION_FORMATS.get(suffix, "")
    if by_ext in {
        "bam",
        "cram",
        "sam",
        "maf",
        "a2m",
        "a3m",
        "selex",
        "mase",
        "xmfa",
        "axt",
        "paf",
        "psl",
        "blast_tab",
        "gaf",
        "hhm",
        "emboss",
        "mega",
        "msf",
        "chain",
        "pdb",
        "mmcif",
        "gbz",
        "odgi",
        "hal",
    }:
        return by_ext
    first = stripped.splitlines()[0].strip() if stripped else ""
    upper = first.upper()
    if upper.startswith("#FORMATVERSION MAUVE"):
        return "xmfa"
    if stripped.startswith(">"):
        return "fasta"
    if upper.startswith("CLUSTAL") or upper.startswith("MUSCLE") or upper.startswith("PROBCONS"):
        return "clustal"
    if upper.startswith("# STOCKHOLM"):
        return "stockholm"
    if upper.startswith("#NEXUS") or stripped.upper().startswith("BEGIN DATA"):
        return "nexus"
    if upper.startswith("#MEGA"):
        return "mega"
    if upper.startswith("HMMER3") or upper.startswith("HMMER2"):
        return "hmmer"
    if upper.startswith("HHSEARCH") or upper.startswith("HH-SUITE"):
        return "hhm"
    if first.startswith("##maf") or (stripped.startswith("a ") and "\ns " in stripped):
        return "maf"
    if (
        "PILEUP" in upper
        or "!!AA_MULTIPLE_ALIGNMENT" in upper
        or "!!NA_MULTIPLE_ALIGNMENT" in upper
    ):
        return "msf"
    if re.search(r"^\s*MSF:\s*\d+", text, re.M) and "//" in text:
        return "msf"
    if first.startswith("########################################") or (
        upper.startswith("#") and "# ALIGNED_SEQUENCES:" in text.upper()
    ):
        return "emboss"
    if upper.startswith("@HD") or upper.startswith("@SQ") or _looks_sam(first):
        return "sam"
    if first[:2] in {"H\t", "S\t", "L\t", "P\t", "W\t"} and _looks_gfa(stripped):
        return "gfa"
    if _looks_paf(first):
        return "paf"
    if suffix in {"phy", "phylip"} or _looks_phylip(first):
        return "phylip"
    if by_ext:
        return by_ext
    raise AlignmentFormatError("Could not auto-detect alignment format.")


def _detect_binary_format(path: Path) -> str:
    suffix = canonical_suffix(path)
    return suffix if suffix in {"bam", "cram"} else "auto"


def _looks_phylip(first_line: str) -> bool:
    parts = first_line.split()
    return len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit()


def _looks_gfa(text: str) -> bool:
    tags = {line.split("\t", 1)[0] for line in text.splitlines()[:50] if "\t" in line}
    return bool(tags) and tags <= {"H", "S", "L", "P", "W", "C", "J", "E", "G", "O", "U", "F"}


def _looks_paf(first_line: str) -> bool:
    parts = first_line.split("\t")
    return (
        len(parts) >= 12
        and all(parts[i].isdigit() for i in (1, 2, 3, 6, 7, 8, 9, 10, 11))
        and parts[4] in {"+", "-"}
    )


def _looks_sam(first_line: str) -> bool:
    parts = first_line.split("\t")
    return len(parts) >= 11 and parts[1].isdigit() and parts[3].isdigit()


def _read_fasta(text: str, keep_case: bool = False) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    current_id: str | None = None
    description = ""
    chunks: list[str] = []

    def flush() -> None:
        seq = "".join(chunks).replace(" ", "").replace("\t", "")
        records.append(
            SequenceRecord(str(current_id), seq if keep_case else seq.upper(), description)
        )

    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if current_id is not None:
                flush()
            header = line[1:].strip()
            current_id = header.split()[0] if header else f"seq_{len(records) + 1}"
            description = header
            chunks = []
        else:
            chunks.append(line)
    if current_id is not None:
        flush()
    if not records:
        raise AlignmentFormatError("No FASTA records found.")
    return records


def _read_phylip(text: str) -> list[SequenceRecord]:
    lines = [line.rstrip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise AlignmentFormatError("Empty PHYLIP file.")
    header = lines[0].split()
    if len(header) < 2 or not header[0].isdigit():
        raise AlignmentFormatError("Invalid PHYLIP header.")
    expected = int(header[0])
    records: list[SequenceRecord] = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 2:
            continue
        seq_id = parts[0]
        seq = "".join(parts[1:]).upper()
        records.append(SequenceRecord(seq_id, seq, seq_id))
        if len(records) == expected:
            break
    if len(records) != expected:
        raise AlignmentFormatError(f"PHYLIP expected {expected} sequences, found {len(records)}.")
    return records


def _read_clustal(text: str) -> list[SequenceRecord]:
    chunks: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.rstrip()
        stripped = line.strip()
        if (
            not stripped
            or stripped.upper().startswith("CLUSTAL")
            or stripped.startswith(("*", ":", "."))
        ):
            continue
        parts = stripped.split()
        if len(parts) >= 2 and not set(parts[0]) <= set("*:."):
            chunks.setdefault(parts[0], []).append(parts[1])
    return [
        SequenceRecord(seq_id, "".join(parts).upper(), seq_id) for seq_id, parts in chunks.items()
    ]


def _read_stockholm(text: str) -> list[SequenceRecord]:
    chunks: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line == "//":
            continue
        parts = line.split()
        if len(parts) >= 2:
            chunks.setdefault(parts[0], []).append(parts[1])
    return [
        SequenceRecord(seq_id, "".join(parts).upper(), seq_id) for seq_id, parts in chunks.items()
    ]


def _read_nexus(text: str) -> list[SequenceRecord]:
    in_matrix = False
    records: list[SequenceRecord] = []
    for raw in text.splitlines():
        line = raw.strip().rstrip(";")
        if not line:
            continue
        if line.upper().startswith("MATRIX"):
            in_matrix = True
            remainder = line[6:].strip()
            if not remainder:
                continue
            line = remainder
        if in_matrix:
            if line.upper().startswith("END") or line == "":
                break
            parts = line.split()
            if len(parts) >= 2 and not parts[0].startswith("["):
                records.append(SequenceRecord(parts[0], "".join(parts[1:]).upper(), parts[0]))
    if not records:
        raise AlignmentFormatError("No NEXUS matrix records found.")
    return records


def _read_maf(text: str) -> list[SequenceRecord]:
    chunks: dict[str, list[str]] = {}
    block_ids: list[str] = []
    n_blocks = 0
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("a"):
            n_blocks += 1
            if n_blocks > 1:
                raise AlignmentFormatError(
                    "This MAF file contains more than one alignment block. MSATK does not "
                    "concatenate blocks into one matrix because that discards genomic "
                    "coordinates; use `msatk inspect`/`msatk profile` (block analysis) or "
                    "msatk.alignments.load_alignment()."
                )
            block_ids = []
            continue
        if line.startswith("s "):
            parts = line.split()
            if len(parts) < 7:
                raise AlignmentFormatError("Invalid MAF sequence line.")
            seq_id = parts[1]
            seq = parts[6].upper()
            block_ids.append(seq_id)
            chunks.setdefault(seq_id, []).append(seq)
    if not chunks:
        raise AlignmentFormatError("No MAF alignment records found.")
    return [SequenceRecord(seq_id, "".join(parts), seq_id) for seq_id, parts in chunks.items()]


def _read_sam(text: str) -> list[SequenceRecord]:
    rows: list[tuple[str, int, str]] = []
    for raw in text.splitlines():
        if not raw.strip() or raw.startswith("@"):
            continue
        parts = raw.rstrip().split("\t")
        if len(parts) < 11:
            continue
        read_id = parts[0]
        pos = int(parts[3])
        cigar = parts[5]
        seq = parts[9]
        if cigar == "*" or seq == "*":
            continue
        rows.append((read_id, pos, _sequence_from_cigar(seq, cigar)))
    return _records_from_reference_offsets(rows, source_name="SAM")


def _read_bam_cram(path: Path, fmt: str) -> list[SequenceRecord]:
    try:
        import pysam  # type: ignore
    except Exception as exc:
        raise AlignmentFormatError(
            "BAM/CRAM input requires the optional 'pysam' dependency. "
            "Install it with conda (`mamba install -c bioconda pysam`) or pip (`pip install pysam`)."
        ) from exc

    mode = "rc" if fmt == "cram" else "rb"
    rows: list[tuple[str, int, str]] = []
    try:
        with pysam.AlignmentFile(str(path), mode) as handle:
            for read in handle.fetch(until_eof=True):
                if read.is_unmapped or read.query_sequence is None or read.cigartuples is None:
                    continue
                seq = _sequence_from_pysam_cigar(read.query_sequence, read.cigartuples)
                rows.append((read.query_name, int(read.reference_start) + 1, seq))
    except Exception as exc:
        raise AlignmentFormatError(
            f"MSATK could not read {fmt.upper()} file {path}: {exc}"
        ) from exc
    return _records_from_reference_offsets(rows, source_name=fmt.upper())


def _records_from_reference_offsets(
    rows: list[tuple[str, int, str]], source_name: str
) -> list[SequenceRecord]:
    if not rows:
        raise AlignmentFormatError(
            f"No mapped reads or alignment records found in {source_name} input."
        )
    min_pos = min(pos for _, pos, _ in rows)
    records = []
    seen: dict[str, int] = {}
    for read_id, pos, seq in rows:
        seen[read_id] = seen.get(read_id, 0) + 1
        unique_id = read_id if seen[read_id] == 1 else f"{read_id}_{seen[read_id]}"
        records.append(SequenceRecord(unique_id, "-" * (pos - min_pos) + seq.upper(), unique_id))
    return records


def _sequence_from_cigar(sequence: str, cigar: str) -> str:
    import re

    pieces: list[str] = []
    query_index = 0
    for length_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_text)
        if op in {"M", "=", "X", "I"}:
            pieces.append(sequence[query_index : query_index + length])
            query_index += length
        elif op in {"D", "N"}:
            pieces.append("-" * length)
        elif op == "S":
            query_index += length
        elif op in {"H", "P"}:
            continue
    return "".join(pieces)


def _sequence_from_pysam_cigar(sequence: str, cigartuples: list[tuple[int, int]]) -> str:
    pieces: list[str] = []
    query_index = 0
    for op, length in cigartuples:
        if op in {0, 1, 7, 8}:  # M, I, =, X
            pieces.append(sequence[query_index : query_index + length])
            query_index += length
        elif op in {2, 3}:  # deletion or skipped region
            pieces.append("-" * length)
        elif op == 4:  # soft clip
            query_index += length
        elif op in {5, 6}:  # hard clip or padding
            continue
    return "".join(pieces)


# --------------------------------------------------------------------------- A2M / A3M
def _read_a2m_a3m(text: str, fmt: str) -> tuple[list[SequenceRecord], list[list[int]]]:
    """Match-state matrix plus per-sequence insert lengths.

    Uppercase letters and ``-`` are match states (aligned columns); lowercase letters
    are insertions relative to the match states. A2M also writes ``.`` for insert
    gaps; A3M omits them. Insert residues are removed from the matrix and counted
    in ``insertions[i][k]`` (residues inserted after match column ``k``; index 0 is
    before the first match column).
    """

    raw = _read_fasta(text, keep_case=True)
    records: list[SequenceRecord] = []
    insertions: list[list[int]] = []
    for rec in raw:
        match: list[str] = []
        ins = [0]
        for ch in rec.sequence:
            if ch.islower():
                ins[-1] += 1
            elif ch == ".":
                continue
            else:
                match.append(ch)
                ins.append(0)
        records.append(SequenceRecord(rec.id, "".join(match).upper(), rec.description))
        insertions.append(ins)
    widths = {len(r.sequence) for r in records}
    if len(widths) > 1:
        raise AlignmentFormatError(
            f"{fmt.upper()} match-state rows have unequal lengths {sorted(widths)}; "
            "the file is not a valid A2M/A3M alignment."
        )
    return records, insertions


# --------------------------------------------------------------------------- Stockholm
def _read_stockholm_annotated(text: str) -> tuple[list[SequenceRecord], dict[str, str]]:
    """Stockholm sequences plus concatenated ``#=GC`` tracks (interleaved blocks)."""

    chunks: dict[str, list[str]] = {}
    gc: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line == "//":
            continue
        if line.startswith("#=GC"):
            parts = line.split(None, 2)
            if len(parts) == 3:
                gc.setdefault(parts[1], []).append(parts[2].replace(" ", ""))
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            chunks.setdefault(parts[0], []).append(parts[1])
    if not chunks:
        raise AlignmentFormatError("No Stockholm sequence records found.")
    records = [
        SequenceRecord(seq_id, "".join(parts).upper(), seq_id) for seq_id, parts in chunks.items()
    ]
    return records, {k: "".join(v) for k, v in gc.items()}


# --------------------------------------------------------------------------- partitions
_RANGE = re.compile(r"(\d+)(?:\s*-\s*(\d+|\.))?(?:\s*\\\s*(\d+))?")


def _expand_ranges(spec: str, width: int | None = None) -> list[int]:
    sites: list[int] = []
    for start, end, step in _RANGE.findall(spec):
        a = int(start)
        b = a if not end else (width or a) if end == "." else int(end)
        sites.extend(range(a, b + 1, int(step) if step else 1))
    return sorted(set(sites))


def parse_nexus_charsets(text: str) -> list[Partition]:
    """Parse ``charset name = 1-300 301-600\\3;`` statements from NEXUS ``sets`` blocks."""

    out: list[Partition] = []
    nchar = re.search(r"NCHAR\s*=\s*(\d+)", text, re.I)
    width = int(nchar.group(1)) if nchar else None
    for name, spec in re.findall(r"charset\s+(\S+?)\s*=\s*([^;]+);", text, re.I):
        out.append(Partition(name.strip("'\""), _expand_ranges(spec, width)))
    return out


def parse_partition_file(path: str | Path, width: int | None = None) -> list[Partition]:
    """Read a RAxML/IQ-TREE style partition file (``DNA, gene1 = 1-300, 301-600\\3``)
    or a NEXUS file containing ``charset`` statements."""

    text = read_text(path)
    if "charset" in text.lower():
        return parse_nexus_charsets(text)
    out: list[Partition] = []
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if "=" not in line:
            continue
        left, spec = line.split("=", 1)
        bits = [b.strip() for b in left.split(",")]
        datatype, name = (bits[0], bits[-1]) if len(bits) > 1 else ("", bits[0])
        out.append(Partition(name, _expand_ranges(spec.replace(",", " "), width), datatype))
    if not out:
        raise AlignmentFormatError(f"No partitions found in {path}.")
    return out


# --------------------------------------------------------------------------- MSF / GCG
def _read_msf(text: str) -> list[SequenceRecord]:
    body = text.split("//", 1)[1] if "//" in text else text
    names = re.findall(r"^\s*Name:\s*(\S+)", text.split("//", 1)[0], re.M)
    chunks: dict[str, list[str]] = {n: [] for n in names}
    for raw in body.splitlines():
        parts = raw.split()
        if len(parts) < 2 or (names and parts[0] not in chunks):
            continue
        if parts[0].isdigit():  # coordinate ruler lines
            continue
        chunks.setdefault(parts[0], []).append("".join(parts[1:]))
    records = [
        SequenceRecord(n, "".join(p).replace("~", "-").replace(".", "-").upper(), n)
        for n, p in chunks.items()
        if p
    ]
    if not records:
        raise AlignmentFormatError("No MSF sequence records found.")
    return records


# --------------------------------------------------------------------------- MEGA
def _read_mega(text: str) -> list[SequenceRecord]:
    """MEGA ``#name sequence`` records; ``.`` (identical to the first sequence) is expanded."""

    chunks: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.upper().startswith("#MEGA") or line.startswith("!"):
            continue
        if line.startswith("#"):
            parts = line[1:].split(None, 1)
            if len(parts) == 2:
                chunks.setdefault(parts[0], []).append(parts[1].replace(" ", ""))
    if not chunks:
        raise AlignmentFormatError("No MEGA sequence records found.")
    seqs = {k: "".join(v).upper() for k, v in chunks.items()}
    first = next(iter(seqs.values()))
    records = []
    for name, seq in seqs.items():
        expanded = "".join(
            first[i] if ch == "." and i < len(first) else ch for i, ch in enumerate(seq)
        )
        records.append(SequenceRecord(name, expanded, name))
    return records


# --------------------------------------------------------------------------- SELEX / MASE
def _read_selex(text: str) -> list[SequenceRecord]:
    chunks: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.rstrip()
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            chunks.setdefault(parts[0], []).append("".join(parts[1:]))
    if not chunks:
        raise AlignmentFormatError("No SELEX sequence records found.")
    return [
        SequenceRecord(n, "".join(p).replace("_", "-").replace(".", "-").upper(), n)
        for n, p in chunks.items()
    ]


def _read_mase(text: str) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    name: str | None = None
    seq: list[str] = []
    expect_name = False
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith(";"):
            if name is not None:
                records.append(SequenceRecord(name, "".join(seq).upper(), name))
                name, seq = None, []
            expect_name = True
            continue
        if not line:
            continue
        if expect_name:
            name, expect_name = line.split()[0], False
        else:
            seq.append(line.replace(" ", ""))
    if name is not None:
        records.append(SequenceRecord(name, "".join(seq).upper(), name))
    if not records:
        raise AlignmentFormatError("No MASE sequence records found.")
    return records


# --------------------------------------------------------------------------- EMBOSS pairwise
_EMBOSS_METHOD = {"needle": "global", "stretcher": "global", "water": "local", "matcher": "local"}


def _read_emboss(text: str) -> tuple[list[SequenceRecord], dict[str, Any]]:
    """First alignment of an EMBOSS pairwise file (needle, water, stretcher, matcher)."""

    meta: dict[str, Any] = {}
    names: list[str] = []
    for raw in text.splitlines():
        m = re.match(r"#\s*(\d):\s*(\S+)", raw)
        if m:
            names.append(m.group(2))
            continue
        m = re.match(r"#\s*(Program|Matrix|Gap_penalty|Extend_penalty|Score|Length):\s*(.+)", raw)
        if m and m.group(1).lower() not in meta:
            meta[m.group(1).lower()] = m.group(2).strip()
    program = str(meta.get("program", "")).lower()
    meta["method"] = _EMBOSS_METHOD.get(program, "unknown")
    rows: list[list[str]] = [[], []]
    starts: list[int | None] = [None, None]
    k = 0
    started = False
    for raw in text.splitlines():
        if raw.startswith("#"):
            if started and raw.startswith("#==="):
                break
            continue
        m = re.match(r"^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s*$", raw)
        if m:
            started = True
            if starts[k] is None:
                starts[k] = int(m.group(2))
            rows[k].append(m.group(3))
            k = 1 - k
    if not rows[0] or len(rows[0]) != len(rows[1]):
        raise AlignmentFormatError("Could not parse an EMBOSS pairwise alignment.")
    names = (names + ["seq1", "seq2"])[:2]
    meta["start_1"], meta["start_2"] = starts
    for key in ("score", "gap_penalty", "extend_penalty"):
        if key in meta:
            with contextlib.suppress(ValueError):
                meta[key] = float(meta[key])
    return [
        SequenceRecord(names[0], "".join(rows[0]).upper(), names[0]),
        SequenceRecord(names[1], "".join(rows[1]).upper(), names[1]),
    ], meta
