"""Characterization analyses: hand-checked values, applicability and registry coverage."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from msatk import MSATK, datasets
from msatk.core.encoded import encode_alignment
from msatk.models import Alignment, SequenceRecord
from msatk.stats.characterization import ANALYSES, PROFILES, characterize, resolve
from msatk.stats.completeness import completeness_analysis, gap_structure_analysis
from msatk.stats.consistency import compare_alignments, ensemble_agreement
from msatk.stats.heterogeneity import _symmetry, composition_heterogeneity
from msatk.stats.pairwise import exact_pairwise
from msatk.stats.redundancy import redundancy_analysis
from msatk.stats.registry import is_defined
from msatk.stats.registry_characterization import TABLE_LEVELS
from msatk.stats.site_patterns import site_pattern_analysis


def aln(*seqs, ids=None):
    ids = ids or [f"s{i}" for i in range(len(seqs))]
    return Alignment(tuple(SequenceRecord(i, s) for i, s in zip(ids, seqs)))


def enc(*seqs, kind="dna"):
    return encode_alignment(aln(*seqs), kind)


# ----------------------------------------------------------------- completeness / gaps
def test_completeness_matches_alistat_definitions():
    e = enc("ACGT", "AC--", "?CGN")
    s, t, m = completeness_analysis(e)
    assert s["alignment_completeness"] == pytest.approx(8 / 12)
    assert s["min_sequence_completeness"] == pytest.approx(0.5)
    assert s["min_site_completeness"] == pytest.approx(1 / 3)  # last column: T only
    # C_ij: (0,1) share A,C -> 2/4; (0,2) share C,G -> 2/4; (1,2) share C -> 1/4
    assert m[0, 1] == pytest.approx(0.5) and m[1, 2] == pytest.approx(0.25)
    assert s["min_pair_completeness"] == pytest.approx(0.25)
    assert s["min_shared_sites"] == 1
    assert s["pairs_below_threshold"] == 1
    row = {r["sequence_id"]: r for r in t["sequence_completeness"]}
    assert row["s1"]["least_overlapping_partner"] == "s2"


def test_gap_structure_counts():
    e = enc("--AC-GT--A", "ACACGGTTAA", "AC---GTTAA")
    s, t = gap_structure_analysis(e)
    rows = {r["sequence_id"]: r for r in t["gap_structure_by_sequence"]}
    # s0: terminal "--" at start; internal runs "-" (len 1) and "--" (len 2)
    assert rows["s0"]["terminal_gap_cells"] == 2
    assert rows["s0"]["internal_gap_openings"] == 2
    assert rows["s0"]["internal_gap_cells"] == 3
    assert rows["s2"]["internal_gap_openings"] == 1 and rows["s2"]["longest_internal_gap"] == 3
    assert s["gap_openings"] == 3 and s["gap_extensions"] == 6 - 3
    assert s["max_gap_run_length"] == 3 and s["single_site_gap_runs"] == 1
    # gap-free columns: 6,7 (G,T) -> one block of 2; also col 10 -> block of 1
    assert s["longest_gap_free_block"] == 2
    dist = {r["run_length"]: r["runs"] for r in t["gap_run_length_distribution"]}
    assert dist == {1: 1, 2: 1, 3: 1}


# ----------------------------------------------------------------- redundancy
def test_neff_identical_and_distinct():
    same = enc("ACGTACGTAC", "ACGTACGTAC", "ACGTACGTAC", "ACGTACGTAT")
    s, t, w = redundancy_analysis(same, exact_pairwise(same))
    assert s["effective_sequences"] == pytest.approx(1.0)  # all four within 80 %
    assert s["clusters"] == 1 and w.sum() == pytest.approx(1.0)
    distinct = enc("AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG")
    s, _, _ = redundancy_analysis(distinct, exact_pairwise(distinct))
    assert s["effective_sequences"] == pytest.approx(3.0)
    assert s["effective_sequence_fraction"] == pytest.approx(1.0)


def test_neff_needs_overlap():
    # Fragment identical to s0 over only 2 of its 10 residues: not a neighbour.
    e = enc("ACGTACGTAC", "AC--------", "TTTTTTTTTT")
    s, t, _ = redundancy_analysis(e, exact_pairwise(e))
    rows = {r["sequence_id"]: r for r in t["sequence_weights"]}
    assert rows["s1"]["neighbours"] == 1
    s2, _, _ = redundancy_analysis(e, exact_pairwise(e), min_overlap=0.0)
    assert s2["effective_sequences"] < s["effective_sequences"]


# ----------------------------------------------------------------- site patterns
def test_site_patterns_and_simpson():
    e = enc("AAAC", "AAAC", "AAGT", "ACGT")
    s, t = site_pattern_analysis(e)
    # columns: AAAA, AAAC, AAGG, CCTT -> 4 distinct
    assert s["distinct_site_patterns"] == 4 and s["sites_1_state"] == 1
    rows = t["site_diversity"]
    assert rows[0]["gini_simpson"] == pytest.approx(0.0)
    assert rows[2]["gini_simpson"] == pytest.approx(0.5)
    assert rows[1]["minor_state_count"] == 1 and rows[2]["minor_state_count"] == 2
    assert s["singleton_minor_sites"] == 1 and s["doubleton_minor_sites"] == 2
    e2 = enc("AAAA", "AAAA", "CCCC")
    s2, t2 = site_pattern_analysis(e2)
    assert s2["distinct_site_patterns"] == 1 and s2["max_pattern_multiplicity"] == 4


# ----------------------------------------------------------------- heterogeneity
def test_rcv_hand_calculation():
    e = enc("AAAA", "CCCC")
    s, t = composition_heterogeneity(e)
    # counts A: 4,0 (mean 2); C: 0,4 (mean 2) -> sum |dev| = 8; RCV = 8 / (2 * 4) = 1
    assert s["rcv"] == pytest.approx(1.0)
    assert sum(r["rcvt"] for r in t["sequence_composition_tests"]) == pytest.approx(1.0)
    homogeneous = enc("ACGT" * 5, "CAGT" * 5, "ACTG" * 5)
    s, t = composition_heterogeneity(homogeneous)
    assert s["rcv"] == pytest.approx(0.0)
    assert all(r["composition_p"] == pytest.approx(1.0) for r in t["sequence_composition_tests"])


def test_symmetry_tests():
    sym = np.array([[10, 3, 1, 0], [3, 12, 2, 1], [1, 2, 9, 4], [0, 1, 4, 11]], float)
    r = _symmetry(sym)
    assert r["bowker"] == pytest.approx(0.0) and r["bowker_p"] == pytest.approx(1.0)
    assert r["marginal"] == pytest.approx(0.0, abs=1e-9)
    asym = sym.copy()
    asym[0, 1], asym[1, 0] = 20, 0
    r = _symmetry(asym)
    assert r["bowker"] == pytest.approx(20.0)  # (20 - 0)^2 / 20
    assert r["bowker_df"] == 5 and r["bowker_p"] < 0.01


# ----------------------------------------------------------------- windows
def test_window_outlier_and_codon_frames():
    rng = np.random.default_rng(0)
    base = "".join(rng.choice(list("ACGT"), size=120))
    seqs = [base] * 12
    divergent = list(base)
    for k in range(60, 90):
        divergent[k] = {"A": "C", "C": "G", "G": "T", "T": "A"}[base[k]]
    seqs.append("".join(divergent))
    e = encode_alignment(aln(*seqs), "codon")
    from msatk.stats.windows import windowed_analysis

    s, t = windowed_analysis(e, window=30, step=15)
    by_start = {r["window_start"]: r["outlier_sequences"] for r in t["windows"]}
    assert by_start[61] == 1 and by_start[1] == 0
    with pytest.raises(ValueError, match="multiples of 3"):
        windowed_analysis(e, window=31)


# ----------------------------------------------------------------- protein properties
def test_protein_conservative_vs_radical():
    from msatk.protein.properties import protein_property_analysis

    e = enc("LLLL", "LLLL", "ILLL", "DLLL", kind="protein")
    s, t = protein_property_analysis(e)
    site = t["site_properties"][0]
    # consensus L; I is conservative (BLOSUM62 L/I = +2), D is radical (L/D = -4)
    assert site["conservative_substitutions"] == 1 and site["radical_substitutions"] == 1
    assert site["class_changes"] == 1  # D is negative, L hydrophobic; I stays hydrophobic
    assert s["conservative_fraction"] == pytest.approx(0.5)


# ----------------------------------------------------------------- codon evolution
def test_ng86_single_changes_and_enc():
    from msatk.codon.evolution import codon_matrix, effective_number_of_codons, pairwise_dnds
    from msatk.core.genetic_code import genetic_code

    a = "CTT" * 12  # Leu
    b = "CTC" + "CTT" * 11  # synonymous change at codon 1
    c = "ATT" + "CTT" * 11  # Leu -> Ile, non-synonymous
    summary, rows = pairwise_dnds(codon_matrix([a, b, c]), ["a", "b", "c"])
    ab = next(r for r in rows if r["sequence_b"] == "b")
    ac = next(r for r in rows if r["sequence_b"] == "c")
    assert ab["synonymous_differences"] == 1 and ab["nonsynonymous_differences"] == 0
    assert ac["synonymous_differences"] == 0 and ac["nonsynonymous_differences"] == 1
    code = genetic_code(1)
    uniform = codon_matrix(["".join(c for c in code.codons if code.codons[c] != "*")] * 4)
    assert effective_number_of_codons(uniform, code) == pytest.approx(61.0)


def test_ng86_matches_biopython_when_stop_paths_allowed():
    pytest.importorskip("Bio.codonalign")
    import warnings

    import msatk.codon.evolution as ev

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
    src = Path(ev.__file__).read_text()
    variant = src.replace(
        """                if code.translate(nxt) == "*":
                    ok = False
                    break
                if code.translate(nxt) == code.translate(cur):""",
        """                if code.translate(nxt) == code.translate(cur):""",
    )
    assert variant != src
    ns: dict = {}
    exec(compile(variant, "ng86_variant", "exec"), ns)
    rng = np.random.default_rng(3)
    stops = {"TAA", "TAG", "TGA"}
    sense = [c for c in ev.CODONS if c not in stops]
    a = [rng.choice(sense) for _ in range(80)]
    b = [x if rng.random() < 0.7 else rng.choice(sense) for x in a]
    _, rows = ns["pairwise_dnds"](ns["codon_matrix"](["".join(a), "".join(b)]), ["a", "b"])
    dn, ds = cal_dn_ds(CodonSeq("".join(a)), CodonSeq("".join(b)), method="NG86")
    assert rows[0]["dn"] == pytest.approx(dn, abs=1e-10)
    assert rows[0]["ds"] == pytest.approx(ds, abs=1e-10)


# ----------------------------------------------------------------- coevolution
def test_mi_apc_recovers_trna_base_pairs():
    from msatk.stats.rna_structure import base_pair_rows

    m = MSATK(datasets.path("trna_rf00005"))
    pairs = {(r["pair_i"], r["pair_j"]) for r in base_pair_rows(m.alignment, m.annotations.ss_cons)}
    top = m.characterize(analyses=["coevolution"]).tables["coevolution_pairs"][:10]
    assert sum((p["site_i"], p["site_j"]) in pairs for p in top) >= 9


# ----------------------------------------------------------------- consistency
def test_compare_identical_and_shifted():
    ref = aln("ACGT-A", "AC-TTA", "ACGTTA")
    s, t = compare_alignments(ref, ref)
    assert s["sp_score"] == 1.0 and s["tc_score"] == 1.0 and s["precision"] == 1.0
    shifted = aln("ACG-TA", "AC-TTA", "ACGTTA")  # s0's T moved one column right
    s, t = compare_alignments(shifted, ref)
    assert 0 < s["sp_score"] < 1 and s["tc_score"] < 1
    bad = {r["site_index"]: r["agreement"] for r in t["column_agreement"]}
    assert bad[5] < 1.0 and bad[1] == 1.0
    e, et = ensemble_agreement(ref, [ref, shifted])
    assert e["alternatives"] == 2 and 0 < e["mean_column_confidence"] < 1
    with pytest.raises(ValueError, match="different residues"):
        compare_alignments(aln("ACGT-A", "AC-TTA", "ACGTTT"), ref)
    with pytest.raises(ValueError, match="different sequences"):
        compare_alignments(aln("ACGT", "ACGT", ids=["x", "y"]), aln("ACGT", "ACGT", ids=["x", "z"]))


# ----------------------------------------------------------------- orchestration
def test_profiles_and_applicability():
    assert set(PROFILES["core"]) <= set(PROFILES["comprehensive"]) <= set(ANALYSES)
    with pytest.raises(ValueError):
        resolve("bogus", None)
    with pytest.raises(ValueError):
        resolve(None, ["nope"])
    m = MSATK(datasets.path("ras_pf00071"))
    r = m.characterize("evolutionary")
    assert r.applicability["codon_evolution"].startswith("skipped")
    assert r.applicability["protein_properties"].startswith("run")
    assert r.applicability["coevolution"].startswith("not requested")
    small = MSATK(datasets.path("ras_pf00071"), max_exact_sequences=10)
    r = small.characterize("core")
    assert "max_exact_sequences" in r.applicability["redundancy"]


@pytest.mark.parametrize("name", ["ras_pf00071", "trna_rf00005", "toy_partitioned"])
def test_every_emitted_field_is_defined(name):
    m = MSATK(datasets.path(name))
    r = m.characterize("all")
    missing = [
        (lvl, k) for lvl, vals in r.summaries.items() for k in vals if not is_defined(k, lvl)
    ]
    for table, rows in r.tables.items():
        level = TABLE_LEVELS.get(table, table)
        missing += [(level, k) for row in rows for k in row if not is_defined(k, level)]
    assert not sorted(set(missing))


def test_codon_fields_defined():
    seqs = [
        "ATGAAACCCGGGTTTAAACCCGGGTTTAAA" + "CTT" * 10,
        "ATGAAGCCAGGGTTCAAACCTGGGTTTAAG" + "CTC" * 10,
        "ATGAAACCCGGATTTAAACCCGGCTTTAAA" + "CTA" * 10,
    ] * 8
    a = aln(*seqs)
    from msatk.stats import compute_statistics

    r = characterize(compute_statistics(a, "codon"), "codon", profile="all")
    assert r.applicability["codon_evolution"].startswith("run")
    missing = [
        (lvl, k) for lvl, vals in r.summaries.items() for k in vals if not is_defined(k, lvl)
    ]
    for table, rows in r.tables.items():
        level = TABLE_LEVELS.get(table, table)
        missing += [(level, k) for row in rows for k in row if not is_defined(k, level)]
    assert not sorted(set(missing))


def test_compute_api():
    m = MSATK(datasets.path("ras_pf00071"))
    out = m.compute(
        [
            "number_of_sequences",
            "redundancy.effective_sequences",
            "rcv",
            "sequence_weights",
            "gap_structure",
        ]
    )
    assert out["number_of_sequences"] == 60
    assert out["redundancy.effective_sequences"] == pytest.approx(60.0)
    assert out["rcv"] > 0
    assert len(out["sequence_weights"]) == 60
    assert "gap_openings" in out["gap_structure"]
    with pytest.raises(KeyError, match="several levels"):
        m.compute(["effective_sequences"])
    with pytest.raises(KeyError, match="Unknown metric"):
        m.compute(["no_such_metric"])
    per_seq = m.compute(["gap_pattern_robust_z"])["gap_pattern_robust_z"]
    assert len(per_seq) == 60 and "sequence_id" in per_seq[0]


def test_profile_writes_characterization(tmp_path):
    m = MSATK(datasets.path("ras_pf00071"))
    m.write_outputs(tmp_path / "out", plots=False, metrics="comprehensive")
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    assert summary["characterization_profile"] == "comprehensive"
    assert "site_patterns" in summary["characterization"]
    for f in [
        "characterization_summary.csv",
        "site_properties.csv",
        "windows.csv",
        "pair_completeness_matrix.csv",
        "sequence_weights.csv",
    ]:
        assert (tmp_path / "out" / "tables" / f).exists(), f
    html = (tmp_path / "out" / "report.html").read_text()
    assert "Alignment characterization" in html
    md = (tmp_path / "out" / "report.md").read_text()
    assert "## Characterization analyses" in md
    assert not math.isnan(summary["characterization"]["redundancy"]["effective_sequences"])


def test_cli_compare_and_metrics(tmp_path, capsys):
    from msatk.cli import main

    ref = tmp_path / "ref.fasta"
    test = tmp_path / "test.fasta"
    ref.write_text(">a\nACGT-A\n>b\nAC-TTA\n>c\nACGTTA\n")
    test.write_text(">a\nACG-TA\n>b\nAC-TTA\n>c\nACGTTA\n")
    assert main(["compare", str(test), "--reference", str(ref), "--out", str(tmp_path / "c")]) == 0
    out = capsys.readouterr().out
    assert "sp_score" in out and (tmp_path / "c" / "tables" / "column_agreement.csv").exists()
    assert main(["compare", str(test), "--alternatives", str(ref), str(test)]) == 0
    assert "mean_column_confidence" in capsys.readouterr().out
    assert main(["metrics", "describe", "redundancy.effective_sequences"]) == 0
    assert "sum_i 1/n_i" in capsys.readouterr().out
    assert main(["metrics", "list", "--level", "coevolution_pairs"]) == 0
    assert "coevolution_pairs.mi_apc" in capsys.readouterr().out
    assert main(["metrics", "describe", "nope"]) == 2
