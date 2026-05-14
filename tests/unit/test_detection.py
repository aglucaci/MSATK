from msatk.core.detect import detect_molecule_type
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
    for name in ["tiny.phy", "tiny.clustal", "tiny.nexus", "tiny.stockholm"]:
        alignment = read_alignment(data_dir / "formats" / name)
        assert alignment.n_sequences == 3
