"""Generate tiny BAM/CRAM fixtures from synthetic reads.

Run from the repository root after installing pysam:

    python tests/data/formats/create_binary_alignment_fixtures.py
"""

from __future__ import annotations

from pathlib import Path


def main() -> None:
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit("Install pysam before generating BAM/CRAM fixtures.") from exc

    root = Path(__file__).parent
    reference = root / "reference.fasta"
    bam_path = root / "tiny.bam"
    cram_path = root / "tiny.cram"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": 10, "SN": "ref"}]}
    reads = [
        ("read1", 0, [(0, 10)], "ATGCTAGCTA"),
        ("read2", 0, [(0, 10)], "ATGCTAGTTA"),
        ("read3", 0, [(0, 3), (2, 1), (0, 6)], "ATGTAGCTA"),
    ]

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as handle:
        for name, start, cigar, sequence in reads:
            segment = pysam.AlignedSegment()
            segment.query_name = name
            segment.query_sequence = sequence
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = start
            segment.mapping_quality = 60
            segment.cigartuples = cigar
            segment.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
            handle.write(segment)
    pysam.index(str(bam_path))

    pysam.faidx(str(reference))
    with (
        pysam.AlignmentFile(str(bam_path), "rb") as source,
        pysam.AlignmentFile(
            str(cram_path), "wc", header=source.header, reference_filename=str(reference)
        ) as target,
    ):
        for read in source.fetch(until_eof=True):
            target.write(read)
    pysam.index(str(cram_path))


if __name__ == "__main__":
    main()
