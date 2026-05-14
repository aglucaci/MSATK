from pathlib import Path

from msatk.core.detect import alignment_detection_summary, detect_molecule_type
from msatk.io.readers import detect_format, read_alignment
from msatk.models import Alignment, SequenceRecord


def test_detect_fasta(data_dir):
    text = (data_dir / "dna" / "tiny_dna.fasta").read_text()
    assert detect_format(text, data_dir / "dna" / "tiny_dna.fasta") == "fasta"


def test_detect_dna():
    alignment = Alignment((SequenceRecord("a", "ATGCTAGCTA"), SequenceRecord("b", "ATGCTAGTTA")))
    assert detect_molecule_type(alignment) == "dna"


def test_detect_protein():
    alignment = Alignment((SequenceRecord("a", "MKTLLV"), SequenceRecord("b", "MKAILV")))
    assert detect_molecule_type(alignment) == "protein"


def test_read_formats(data_dir):
    for name in [
        "tiny.aligned.fasta",
        "tiny.a3m",
        "tiny.phy",
        "tiny_relaxed.phylip",
        "tiny.clustal",
        "tiny.nexus",
        "tiny.stockholm",
        "tiny.maf",
        "tiny.sam",
    ]:
        alignment = read_alignment(data_dir / "formats" / name)
        assert alignment.n_sequences == 3


def test_detect_binary_alignment_suffixes():
    assert detect_format("", Path("reads.bam")) == "bam"
    assert detect_format("", Path("reads.cram")) == "cram"


def test_alignment_detection_summary_reports_codon_flags(data_dir):
    alignment = read_alignment(data_dir / "codon" / "codon_with_stop.fasta")
    summary = alignment_detection_summary(alignment, "codon")
    assert summary["detected_file_format"] == "fasta"
    assert summary["length_divisible_by_three"] is True
    assert summary["stop_codons_exist"] is True
    assert summary["internal_stop_codons_exist"] is True
    assert summary["sequence_lengths_consistent"] is True
