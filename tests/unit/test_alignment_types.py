"""Alignment families, formats, capability registry and specialized analyses."""

from __future__ import annotations

import gzip
import importlib.util
import math
import shutil
from pathlib import Path

import pytest

from msatk import MSATK
from msatk.alignments import (
    FORMATS,
    BlockAlignment,
    GraphAlignment,
    MappingAlignment,
    ProfileModel,
    ReadAlignment,
    SequenceAlignment,
    compatibility_markdown,
    inspect_alignment,
    load_alignment,
)
from msatk.core.genetic_code import available_codes, genetic_code
from msatk.exceptions import AlignmentFormatError
from msatk.io.readers import parse_partition_file, read_alignment, read_sequence_alignment
from msatk.stats.registry import is_defined

ROOT = Path(__file__).resolve().parents[2]
FMT = ROOT / "tests" / "data" / "formats"
HAS_PYSAM = importlib.util.find_spec("pysam") is not None


# --------------------------------------------------------------------------- registry
def test_registry_fixtures_exist_and_load():
    for cap in FORMATS:
        if cap.level in {"planned", "experimental"}:
            assert not cap.tested_with or cap.level == "experimental"
            continue
        assert cap.tested_with, f"{cap.name} has no test fixture"
        for rel in cap.tested_with:
            path = ROOT / rel
            assert path.exists(), rel
            if cap.dependency == "pysam" and not HAS_PYSAM:
                continue
            obj = load_alignment(path, cap.name)
            assert obj.kind == cap.kind, (rel, obj.kind, cap.kind)


def test_compatibility_page_is_current():
    page = ROOT / "docs" / "compatibility.md"
    assert page.read_text().strip() == compatibility_markdown().strip(), (
        "Regenerate with: msatk formats > docs/compatibility.md"
    )


def test_planned_format_is_reported_not_parsed(tmp_path):
    pdb = tmp_path / "x.pdb"
    pdb.write_text(
        "ATOM      1  N   MET A   1      11.104  13.207   2.100  1.00  0.00           N\n"
    )
    info = inspect_alignment(pdb)
    assert info["kind"] == "structure" and info["support_level"] == "planned"
    assert info["status"].startswith("not readable")
    with pytest.raises(AlignmentFormatError, match="not supported yet"):
        load_alignment(pdb)


@pytest.mark.parametrize("name", ["tiny.aligned.fasta", "blocks.maf", "reads.sam", "tiny.gfa"])
def test_gzip_input(tmp_path, name):
    target = tmp_path / (name + ".gz")
    with (FMT / name).open("rb") as src, gzip.open(target, "wb") as dst:
        shutil.copyfileobj(src, dst)
    plain = load_alignment(FMT / name)
    packed = load_alignment(target)
    assert packed.kind == plain.kind
    assert inspect_alignment(target)["compressed"] == "gzip"


# --------------------------------------------------------------------------- sequence formats
@pytest.mark.parametrize("name", ["tiny.msf", "tiny.meg", "tiny.slx", "tiny.mase"])
def test_specialized_msa_formats_match_fasta(name):
    ref = read_alignment(FMT / "tiny.aligned.fasta")
    other = read_alignment(FMT / name)
    assert other.ids == ref.ids
    assert other.sequences == ref.sequences


@pytest.mark.parametrize("name", ["tiny.a3m", "tiny.a2m"])
def test_a2m_a3m_match_and_insert_states(name):
    aln, ann = read_sequence_alignment(FMT / name)
    assert aln.sequences == ["ATGTAGCTA", "ATGTAGTTA", "ATGTAGCTA"]
    assert ann.insertions is not None
    assert [sum(v) for v in ann.insertions] == [1, 1, 0]
    assert ann.insertions[0][3] == 1  # inserted after the third match state


def test_emboss_pairwise_statistics_match_reported_score():
    obj = load_alignment(FMT / "tiny.needle")
    assert isinstance(obj, SequenceAlignment)
    assert obj.scope == "pairwise" and obj.annotations.method == "global"
    s = obj.summary()["specialized"]["pairwise_alignment"]
    assert s["identical_pairs"] == 6 and s["mismatched_pairs"] == 1
    assert s["identity"] == pytest.approx(6 / 7)
    assert s["identity_over_columns"] == pytest.approx(6 / 9)
    assert s["similarity"] == pytest.approx(1.0)  # I/L is BLOSUM62-positive
    assert s["gap_openings"] == 1 and s["gap_extensions"] == 1
    assert s["score_recomputed"] == pytest.approx(19.5) == s["score_reported"]


def test_multi_block_maf_is_not_concatenated():
    with pytest.raises(AlignmentFormatError, match="more than one alignment block"):
        read_alignment(FMT / "blocks.maf")


# --------------------------------------------------------------------------- blocks
def test_maf_block_statistics():
    obj = load_alignment(FMT / "blocks.maf")
    assert isinstance(obj, BlockAlignment)
    s = obj.summary()
    assert s["blocks"] == 3 and s["genomes"] == 3 and s["reference_genome"] == "hg38"
    assert s["reference_aligned_fraction"] == pytest.approx(24 / 1000)
    assert s["minus_strand_rows"] == 1
    assert s["collinearity_breakpoints"] == 1  # mm10 block 3 is inverted
    genomes = {g["genome"]: g for g in obj.genome_rows()}
    assert genomes["canFam"]["blocks_present"] == 2
    assert genomes["mm10"]["covered_bases"] == 24
    # minus-strand MAF start is converted to forward coordinates: 2000 - 1400 - 6
    assert obj.blocks[2].rows[1].start == 594


def test_xmfa_and_axt_blocks():
    x = load_alignment(FMT / "tiny.xmfa")
    assert isinstance(x, BlockAlignment)
    assert x.genomes == ["genomeA", "genomeB"]  # file names, extension dropped
    assert len(x.blocks[1].rows) == 1  # 2:0-0 marks an absent genome
    assert x.blocks[0].rows[1].strand == "-"
    a = load_alignment(FMT / "tiny.axt")
    assert isinstance(a, BlockAlignment) and len(a.blocks) == 2


# --------------------------------------------------------------------------- reads
def test_sam_read_statistics():
    obj = load_alignment(FMT / "reads.sam")
    assert isinstance(obj, ReadAlignment)
    s = obj.summary()
    assert s["records"] == 8 and s["primary_reads"] == 7 and s["secondary_records"] == 1
    assert s["mapping_rate"] == pytest.approx(6 / 7)
    assert s["mean_mapq"] == pytest.approx(45)
    assert s["aligned_bases"] == 93
    assert s["inserted_bases"] == 2 and s["deleted_bases"] == 5
    assert s["mismatch_rate"] == pytest.approx(2 / 93)  # NM 9 − I 2 − D 5
    assert s["spliced_reads"] == 1 and s["soft_clipped_bases"] == 5
    assert s["median_insert_size"] == 50
    refs = {r["reference"]: r for r in obj.reference_rows()}
    assert refs["chrA"]["breadth_1x"] == pytest.approx(53 / 100)  # duplicate r6 excluded
    assert refs["chrA"]["mean_depth"] == pytest.approx(63 / 100)
    assert refs["chrB"]["breadth_1x"] == pytest.approx(20 / 50)
    assert s["breadth_1x"] == pytest.approx(73 / 150)


@pytest.mark.skipif(not HAS_PYSAM, reason="pysam is required")
def test_bam_read_statistics():
    obj = load_alignment(FMT / "tiny.bam")
    assert obj.kind == "read" and obj.summary()["records"] >= 1


# --------------------------------------------------------------------------- mappings
def test_paf_psl_blast_gaf():
    paf = load_alignment(FMT / "tiny.paf")
    assert isinstance(paf, MappingAlignment)
    s = paf.summary()
    assert s["mean_identity"] == pytest.approx((0.99 + 0.9 + 700 / 750) / 3)
    assert s["multi_hit_queries"] == 1
    assert s["mean_query_coverage"] == pytest.approx((1.0 + 750 / 800) / 2)
    psl = load_alignment(FMT / "tiny.psl").summary()
    assert psl["mean_identity"] == pytest.approx((0.97 + 0.8) / 2)
    blast = load_alignment(FMT / "tiny.m8")
    assert [r.strand for r in blast.records] == ["+", "-", "+"]
    assert blast.summary()["min_evalue"] == pytest.approx(1e-100)
    gaf = load_alignment(FMT / "tiny.gaf")
    assert gaf.fmt == "gaf" and gaf.records[0].target == ">s1>s2>s4"


# --------------------------------------------------------------------------- graphs
def test_gfa_pangenome_statistics():
    obj = load_alignment(FMT / "tiny.gfa")
    assert isinstance(obj, GraphAlignment)
    s = obj.summary()
    assert s["segments"] == 5 and s["links"] == 4 and s["connected_components"] == 2
    assert s["tip_segments"] == 2 and s["isolated_segments"] == 1
    assert (s["core_segments"], s["shell_segments"], s["private_segments"]) == (2, 1, 1)
    assert s["core_bases_fraction"] == pytest.approx(15 / 24)
    assert s["reference_rank0_segments"] == 3


def test_gfa2_rejected(tmp_path):
    p = tmp_path / "g.gfa"
    p.write_text("H\tVN:Z:2.0\nS\ts1\t3\tACG\n")
    with pytest.raises(AlignmentFormatError, match="GFA 2"):
        load_alignment(p)


# --------------------------------------------------------------------------- profiles
def test_hmmer_and_hhm_profiles():
    hmm = load_alignment(FMT / "tiny.hmm")
    assert isinstance(hmm, ProfileModel)
    s = hmm.summary()
    assert s["model_length"] == 3 and s["alphabet"] == "DNA"
    h1 = -(0.85 * math.log2(0.85) + 3 * 0.05 * math.log2(0.05))
    assert s["mean_match_entropy_bits"] == pytest.approx((h1 + 2.0 + 0.0) / 3, abs=1e-4)
    hhm = load_alignment(FMT / "tiny.hhm").summary()
    assert hhm["model_length"] == 4 and hhm["effective_sequences"] == pytest.approx(3.2)


# --------------------------------------------------------------------------- specialized
def test_rna_structure_covariation():
    obj = load_alignment(FMT / "rna_structure.sto")
    special = obj.summary(pairwise="none")["specialized"]
    rna = special["rna_structure"]
    assert rna["consensus_base_pairs"] == 3
    assert rna["base_pairs_with_compensatory_change"] == 2  # GC->CG at 1-9 and 3-7
    assert rna["mean_canonical_pair_fraction"] == pytest.approx(1.0)
    prof = special["profile_states"]
    assert prof["insert_columns"] == 1 and prof["total_inserted_residues"] == 4


def test_partitions_from_nexus_and_file(tmp_path):
    obj = load_alignment(FMT / "partitioned.nexus")
    names = [p.name for p in obj.annotations.partitions]
    assert names == ["geneA", "geneB", "pos3"]
    assert obj.annotations.partitions[2].sites == [3, 6, 9, 12]
    part = tmp_path / "parts.txt"
    part.write_text("DNA, first = 1-6\nDNA, second = 7-12\n")
    parts = parse_partition_file(part, 12)
    assert [p.sites for p in parts] == [list(range(1, 7)), list(range(7, 13))]
    prof = MSATK(FMT / "partitioned.nexus", sequence_type="dna", partitions=part)
    s = prof.specialized()[0]["partitions"]
    assert s["partitions"] == 2 and s["partition_sites_unassigned"] == 0
    rows = prof.specialized()[1]["partitions"]
    assert rows[0]["sequences_without_data"] == 0 and rows[0]["variable_sites"] == 1


def test_composition_chi2_matches_reference_value():
    from msatk.stats.partitions import _gammainc_upper

    # scipy.stats.chi2.sf(3.7916666666666665, 9) = 0.924581591380513
    assert _gammainc_upper(4.5, 3.7916666666666665 / 2) == pytest.approx(0.924581591380513)
    assert _gammainc_upper(1.5, 100) == pytest.approx(4.218541107192018e-43, rel=1e-9)


def test_blosum62_column_scores():
    from msatk.core.encoded import encode_alignment
    from msatk.models import Alignment, SequenceRecord
    from msatk.stats.substitution import blosum62_score, column_substitution_scores

    aln = Alignment(
        (SequenceRecord("a", "WL-"), SequenceRecord("b", "WI-"), SequenceRecord("c", "WLA"))
    )
    rows = column_substitution_scores(encode_alignment(aln, "protein"))
    assert rows[0]["mean_blosum62"] == pytest.approx(11)  # W/W
    assert rows[1]["mean_blosum62"] == pytest.approx((2 + 4 + 2) / 3)  # L/I, L/L, I/L
    assert rows[2]["mean_blosum62"] is None
    assert blosum62_score("W", "W") == 11 and blosum62_score("L", "I") == 2


# --------------------------------------------------------------------------- genetic codes
def test_genetic_codes_bundle_all_ncbi_tables():
    codes = available_codes()
    assert {1, 2, 4, 5, 11, 33} <= set(codes)
    assert genetic_code(1).translate("TGA") == "*"
    assert genetic_code(2).translate("TGA") == "W"  # vertebrate mitochondrial
    assert {"AGA", "AGG"} <= genetic_code(2).stops
    assert genetic_code(4).translate("TGA") == "W"  # mycoplasma
    with pytest.raises(ValueError):
        genetic_code(99)


def test_genetic_code_changes_stop_codon_qc(tmp_path):
    p = tmp_path / "mito.fasta"
    p.write_text(">a\nATGTGAAAATAA\n>b\nATGTGGAAATAA\n")
    std = MSATK(p, sequence_type="codon")
    mito = MSATK(p, sequence_type="codon", genetic_code=2)
    assert any(f.code == "INTERNAL_STOP_CODON" for f in std.qc_findings())
    assert not any(f.code == "INTERNAL_STOP_CODON" for f in mito.qc_findings())
    assert std.stats().summary["genetic_code"] == 1 and mito.stats().summary["genetic_code"] == 2
    from msatk.codon.analysis import rscu

    families = {r["amino_acid"]: r["synonymous_family_size"] for r in rscu(mito.alignment, 2)}
    assert families["W"] == 2 and families["R"] == 4  # AGA/AGG are stops in code 2


# --------------------------------------------------------------------------- metric registry
@pytest.mark.parametrize(
    "name,level,table_levels",
    [
        ("blocks.maf", "block", {"blocks": "block_row", "genomes": "genome"}),
        ("reads.sam", "read", {"references": "reference"}),
        ("tiny.paf", "mapping", {"queries": "mapping_interval", "targets": "mapping_interval"}),
        ("tiny.gfa", "graph", {"segments": "segment", "paths": "path"}),
        ("tiny.hmm", "profile_model", {"match_states": "match_state"}),
    ],
)
def test_kind_outputs_are_registered(name, level, table_levels):
    obj = load_alignment(FMT / name)
    missing = [k for k in obj.summary() if not is_defined(k, level)]
    for table, rows in obj.tables().items():
        if table in table_levels:
            missing += [k for r in rows for k in r if not is_defined(k, table_levels[table])]
    assert not sorted(set(missing))


@pytest.mark.parametrize(
    "name", ["tiny.needle", "rna_structure.sto", "partitioned.nexus", "tiny.a3m"]
)
def test_specialized_outputs_are_registered(name):
    prof = MSATK(FMT / name)
    summaries, tables, _ = prof.specialized()
    levels = {
        "rna_base_pairs": "rna_pair",
        "profile_states": "profile_sequence",
        "partitions": "partition",
        "site_substitution_scores": "site_substitution",
    }
    missing = [(lvl, k) for lvl, vals in summaries.items() for k in vals if not is_defined(k, lvl)]
    missing += [
        (t, k)
        for t, rows in tables.items()
        for r in rows
        for k in r
        if not is_defined(k, levels[t])
    ]
    assert not missing
    summary = prof._provenance_summary()
    assert not [k for k in summary if not is_defined(k, "alignment")]
