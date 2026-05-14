"""Input format detection and alignment readers."""

from __future__ import annotations

from pathlib import Path

from msatk.exceptions import AlignmentFormatError
from msatk.io.validate import validate_alignment
from msatk.models import Alignment, SequenceRecord


def read_alignment(
    path: str | Path, fmt: str = "auto", validation_mode: str = "permissive"
) -> Alignment:
    """Read an alignment from FASTA/A3M, PHYLIP, CLUSTAL, Stockholm, or NEXUS."""

    source = Path(path)
    text = source.read_text(encoding="utf-8")
    if not text.strip():
        raise AlignmentFormatError(f"MSATK could not read {source}: the file is empty.")
    detected = detect_format(text, source) if fmt == "auto" else fmt.lower()
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
    else:
        raise AlignmentFormatError(f"Unsupported alignment format: {fmt}")
    alignment = Alignment(tuple(records), source=str(source), fmt=detected)
    validate_alignment(alignment, mode=validation_mode)
    return alignment


def detect_format(text: str, path: Path | None = None) -> str:
    stripped = text.lstrip()
    suffix = path.suffix.lower().lstrip(".") if path else ""
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
    if suffix in {"phy", "phylip"} or _looks_phylip(first):
        return "phylip"
    if suffix in {"fa", "fasta", "faa", "fna", "a3m", "aln", "sto", "nex"}:
        return {
            "fa": "fasta",
            "faa": "fasta",
            "fna": "fasta",
            "aln": "clustal",
            "sto": "stockholm",
            "nex": "nexus",
        }.get(suffix, suffix)
    raise AlignmentFormatError("Could not auto-detect alignment format.")


def _looks_phylip(first_line: str) -> bool:
    parts = first_line.split()
    return len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit()


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
