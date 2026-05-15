"""Input format detection and alignment readers."""

from __future__ import annotations

from pathlib import Path

from msatk.exceptions import AlignmentFormatError
from msatk.io.validate import validate_alignment
from msatk.models import Alignment, SequenceRecord


def read_alignment(
    path: str | Path, fmt: str = "auto", validation_mode: str = "permissive"
) -> Alignment:
    """Read common MSA formats plus SAM/BAM/CRAM-derived alignments."""

    source = Path(path)
    detected = _detect_binary_format(source) if fmt == "auto" else fmt.lower()
    text = ""
    if detected not in {"bam", "cram"}:
        text = source.read_text(encoding="utf-8")
        if not text.strip():
            raise AlignmentFormatError(f"MSATK could not read {source}: the file is empty.")
        detected = detect_format(text, source) if fmt == "auto" else detected
    if detected in {"fasta", "fa", "faa", "fna", "a3m"}:
        records = _read_fasta(text, keep_lower=detected != "a3m")
    elif detected in {"phylip", "phy"}:
        records = _read_phylip(text)
    elif detected in {"clustal", "aln"}:
        records = _read_clustal(text)
    elif detected in {"stockholm", "sto"}:
        records = _read_stockholm(text)
    elif detected in {"nexus", "nex"}:
        records = _read_nexus(text)
    elif detected == "maf":
        records = _read_maf(text)
    elif detected == "sam":
        records = _read_sam(text)
    elif detected in {"bam", "cram"}:
        records = _read_bam_cram(source, detected)
    else:
        raise AlignmentFormatError(f"Unsupported alignment format: {fmt}")
    alignment = Alignment(tuple(records), source=str(source), fmt=detected)
    validate_alignment(alignment, mode=validation_mode)
    return alignment


def detect_format(text: str, path: Path | None = None) -> str:
    stripped = text.lstrip()
    suffix = path.suffix.lower().lstrip(".") if path else ""
    if suffix in {"bam", "cram", "sam", "maf"}:
        return suffix
    if stripped.startswith(">"):
        return "a3m" if suffix == "a3m" else "fasta"
    first = stripped.splitlines()[0].strip() if stripped else ""
    upper = first.upper()
    if upper.startswith("CLUSTAL"):
        return "clustal"
    if upper.startswith("# STOCKHOLM"):
        return "stockholm"
    if upper.startswith("#NEXUS") or stripped.upper().startswith("BEGIN DATA"):
        return "nexus"
    if stripped.startswith("a ") and "\ns " in stripped:
        return "maf"
    if upper.startswith("@HD") or upper.startswith("@SQ") or _looks_sam(first):
        return "sam"
    if suffix in {"phy", "phylip"} or _looks_phylip(first):
        return "phylip"
    if suffix in {"fa", "fasta", "faa", "fna", "a3m", "aln", "sto", "nex", "maf", "sam"}:
        return {
            "fa": "fasta",
            "faa": "fasta",
            "fna": "fasta",
            "aln": "clustal",
            "sto": "stockholm",
            "nex": "nexus",
        }.get(suffix, suffix)
    raise AlignmentFormatError("Could not auto-detect alignment format.")


def _detect_binary_format(path: Path) -> str:
    suffix = path.suffix.lower().lstrip(".")
    return suffix if suffix in {"bam", "cram"} else "auto"


def _looks_phylip(first_line: str) -> bool:
    parts = first_line.split()
    return len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit()


def _looks_sam(first_line: str) -> bool:
    parts = first_line.split("\t")
    return len(parts) >= 11 and parts[1].isdigit() and parts[3].isdigit()


def _read_fasta(text: str, keep_lower: bool = True) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    current_id: str | None = None
    description = ""
    chunks: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records.append(
                    SequenceRecord(current_id, _clean_sequence(chunks, keep_lower), description)
                )
            header = line[1:].strip()
            current_id = header.split()[0] if header else f"seq_{len(records) + 1}"
            description = header
            chunks = []
        else:
            chunks.append(line)
    if current_id is not None:
        records.append(SequenceRecord(current_id, _clean_sequence(chunks, keep_lower), description))
    if not records:
        raise AlignmentFormatError("No FASTA records found.")
    return records


def _clean_sequence(chunks: list[str], keep_lower: bool = True) -> str:
    sequence = "".join(chunks).replace(" ", "").replace("\t", "")
    if not keep_lower:
        sequence = "".join(ch for ch in sequence if not ch.islower())
    return sequence.upper()


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
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("a "):
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
