"""Generate a tiny BAM fixture without external dependencies."""

from __future__ import annotations

import binascii
import struct
import zlib
from pathlib import Path

BASE_CODES = {
    "=": 0,
    "A": 1,
    "C": 2,
    "M": 3,
    "G": 4,
    "R": 5,
    "S": 6,
    "V": 7,
    "T": 8,
    "W": 9,
    "Y": 10,
    "H": 11,
    "K": 12,
    "D": 13,
    "B": 14,
    "N": 15,
}


def main() -> None:
    root = Path(__file__).parent
    header_text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:10\n"
    payload = bytearray()
    payload.extend(b"BAM\1")
    payload.extend(struct.pack("<i", len(header_text)))
    payload.extend(header_text.encode())
    payload.extend(struct.pack("<i", 1))
    payload.extend(struct.pack("<i", 4))
    payload.extend(b"ref\0")
    payload.extend(struct.pack("<i", 10))
    reads = [
        ("read1", 0, [(10, 0)], "ATGCTAGCTA"),
        ("read2", 0, [(10, 0)], "ATGCTAGTTA"),
        ("read3", 0, [(3, 0), (1, 2), (6, 0)], "ATGTAGCTA"),
    ]
    for read in reads:
        payload.extend(_alignment_record(*read))
    (root / "tiny.bam").write_bytes(_bgzf_block(bytes(payload)) + BGZF_EOF)


def _alignment_record(name: str, pos: int, cigar: list[tuple[int, int]], seq: str) -> bytes:
    read_name = f"{name}\0".encode()
    encoded_seq = _encode_sequence(seq)
    qual = bytes([30] * len(seq))
    n_cigar = len(cigar)
    l_read_name = len(read_name)
    bin_mq_nl = (4681 << 16) | (60 << 8) | l_read_name
    flag_nc = n_cigar
    body = bytearray()
    body.extend(struct.pack("<iiIIiiii", 0, pos, bin_mq_nl, flag_nc, len(seq), -1, -1, 0))
    body.extend(read_name)
    for length, op in cigar:
        body.extend(struct.pack("<I", (length << 4) | op))
    body.extend(encoded_seq)
    body.extend(qual)
    return struct.pack("<i", len(body)) + bytes(body)


def _encode_sequence(seq: str) -> bytes:
    values = [BASE_CODES.get(base.upper(), 15) for base in seq]
    packed = bytearray()
    for index in range(0, len(values), 2):
        high = values[index]
        low = values[index + 1] if index + 1 < len(values) else 0
        packed.append((high << 4) | low)
    return bytes(packed)


def _bgzf_block(data: bytes) -> bytes:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(data) + compressor.flush()
    total_size = 18 + len(compressed) + 8
    header = (
        b"\x1f\x8b\x08\x04"
        + struct.pack("<I", 0)
        + b"\x00\xff"
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<HH", 2, total_size - 1)
    )
    footer = struct.pack("<II", binascii.crc32(data) & 0xFFFFFFFF, len(data) & 0xFFFFFFFF)
    return header + compressed + footer


BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


if __name__ == "__main__":
    main()
