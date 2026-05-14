def test_supported_text_and_binary_fixture_inventory(data_dir):
    formats = data_dir / "formats"
    expected = [
        "tiny.aligned.fasta",
        "tiny.a3m",
        "tiny.phy",
        "tiny_relaxed.phylip",
        "tiny.clustal",
        "tiny.stockholm",
        "tiny.nexus",
        "tiny.maf",
        "tiny.sam",
        "tiny.bam",
        "reference.fasta",
        "create_binary_alignment_fixtures.py",
    ]
    for name in expected:
        assert (formats / name).exists(), name
