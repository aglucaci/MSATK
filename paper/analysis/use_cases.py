"""Reproduce the two Application Note use cases and the figure panels they feed.

(i)  Influenza A(H3N2) haemagglutinin coding sequences (400 x 1,701 nt, codon):
     frame integrity, codon-position entropy, QC findings cross-tabulated with
     host labels from the GenBank titles, and the effect of QC filtering on
     population-genetic summaries.
(ii) Pfam PF00071 (Ras family) seed alignment (60 x 229 aa): conservation of the
     G-box motifs, physicochemical variability, residue-level figures and a
     windowed one-hot export loaded through a PyTorch ``DataLoader``.

Usage::

    python paper/analysis/prepare_data.py --offline
    python paper/analysis/use_cases.py

Writes ``paper/analysis/results/use_cases.json`` and panel images in
``paper/analysis/results/panels/``.
"""

from __future__ import annotations

import contextlib
import csv
import json
import re
import statistics
import time
from pathlib import Path

import numpy as np

from msatk import MSATK
from msatk.encoding.tensors import categorical, decode_categorical, physicochemical
from msatk.encoding.transforms import subset_alignment
from msatk.encoding.vocab import full_vocabulary
from msatk.ml.datasets import AlignmentTensorDataset, to_torch
from msatk.stats import compute_statistics

HERE = Path(__file__).resolve().parent
with contextlib.suppress(ValueError):  # record repo-relative paths in manifests
    HERE = HERE.relative_to(Path.cwd())
DATA = HERE / "data"
RESULTS = HERE / "results"
PANELS = RESULTS / "panels"

G_BOXES = {
    "G1 (P-loop)": r"G..G.GK[ST]",
    "G3 (DxxG)": r"D[TS]AG[QL]E",
    "G4 (NKxD)": r"[NT]K.D",
    "G5 (SAK)": r"[ETS]A[KR]",
}


def _r(x: float | None, nd: int = 4) -> float | None:
    return None if x is None else round(float(x), nd)


def influenza() -> dict[str, object]:
    t0 = time.perf_counter()
    profiler = MSATK(DATA / "h3n2_ha.fasta", sequence_type="codon")
    profiler.write_outputs(RESULTS / "h3n2", preset="comprehensive", force=True)
    seconds = time.perf_counter() - t0
    s = json.loads((RESULTS / "h3n2" / "summary.json").read_text())
    findings = json.loads((RESULTS / "h3n2" / "qc_findings.json").read_text())
    findings = findings if isinstance(findings, list) else findings["findings"]
    flagged = {f["code"]: set(f.get("items") or []) for f in findings}
    with (DATA / "h3n2_ha_metadata.tsv").open() as handle:
        host = {r["accession"]: r["host_label"] for r in csv.DictReader(handle, delimiter="\t")}
    per_seq = {r["sequence_id"]: r for r in profiler.stats().sequences}
    nonhuman = {a for a, h in host.items() if h != "human"}
    divergent = flagged.get("DIVERGENT_SEQUENCES", set())
    unflagged_nonhuman = sorted(nonhuman - divergent)
    with (RESULTS / "h3n2" / "tables" / "codon_position_entropy.csv").open() as handle:
        cpe = {int(r["codon_position"]): float(r["mean_entropy"]) for r in csv.DictReader(handle)}

    # Population-genetic summaries before and after removing QC-flagged divergent sequences.
    keep = [a for a in profiler.alignment.ids if a not in divergent]
    filtered, _ = subset_alignment(profiler.alignment, ids=keep)
    fs = compute_statistics(filtered, "codon", pairwise="exact").summary
    popgen_keys = [
        "nucleotide_diversity_pi",
        "theta_watterson",
        "tajimas_d",
        "haplotype_diversity",
        "ts_tv_ratio",
        "variable_sites",
        "mean_pairwise_identity",
    ]
    return {
        "dataset": "Influenza A(H3N2) HA CDS, NCBI nuccore (2017 publication date), 400 records",
        "runtime_seconds_comprehensive_profile": round(seconds, 2),
        "sequences": s["number_of_sequences"],
        "sites": s["alignment_length"],
        "codons": s["alignment_length"] // 3,
        "unique_sequences": s["unique_sequence_count"],
        "duplicate_sequences": s["duplicate_sequence_count"],
        "terminal_stop_codons": s["stop_codon_count"],
        "internal_stop_codons": s["internal_stop_codons_exist"],
        "frameshift_warnings": s["frameshift_warning_count"],
        "ambiguous_cells": s["ambiguous_cells"],
        "codon_position_mean_entropy": cpe,
        "qc_findings": [
            {"code": f["code"], "severity": f["severity"], "count": f.get("count")}
            for f in findings
        ],
        "host_labels": {
            h: sum(1 for v in host.values() if v == h) for h in sorted(set(host.values()))
        },
        "nonhuman_sequences": len(nonhuman),
        "nonhuman_flagged_divergent": len(nonhuman & divergent),
        "human_flagged_divergent": len(divergent - nonhuman),
        "unflagged_nonhuman": [
            {
                "accession": a,
                "host": host[a],
                "consensus_identity": _r(per_seq[a]["consensus_identity"]),
                "nearest_neighbor_host": host.get(str(per_seq[a]["nearest_neighbor"])),
                "nearest_neighbor_identity": _r(per_seq[a]["nearest_neighbor_identity"]),
            }
            for a in unflagged_nonhuman
        ],
        "consensus_identity_range": {
            "flagged_nonhuman": [
                _r(min(float(per_seq[a]["consensus_identity"]) for a in nonhuman & divergent)),
                _r(max(float(per_seq[a]["consensus_identity"]) for a in nonhuman & divergent)),
            ],
            "unflagged": [
                _r(min(float(per_seq[a]["consensus_identity"]) for a in keep)),
                _r(max(float(per_seq[a]["consensus_identity"]) for a in keep)),
            ],
        },
        "popgen_all": {k: _r(s.get(k)) for k in popgen_keys},
        "popgen_after_qc": {k: _r(fs.get(k)) for k in popgen_keys} | {"sequences": len(keep)},
    }


def ras() -> dict[str, object]:
    path = DATA / "PF00071.fasta"
    t0 = time.perf_counter()
    profiler = MSATK(path, sequence_type="protein", pairwise="exact")
    bundle = profiler.stats()
    stats_seconds = time.perf_counter() - t0
    s, sites = bundle.summary, bundle.sites
    ref = next(r for r in profiler.alignment.records if r.id.startswith("RAS2_YEAST"))
    cols = [i for i, c in enumerate(ref.sequence) if c != "-"]
    ungapped = "".join(c for c in ref.sequence if c != "-")
    motifs = {}
    motif_sites: set[int] = set()
    for name, pattern in G_BOXES.items():
        m = re.search(pattern, ungapped)
        if not m:
            continue
        span = list(range(cols[m.start()] + 1, cols[m.end() - 1] + 2))
        motif_sites.update(span)
        motifs[name] = {
            "reference_motif": m.group(),
            "sites": [span[0], span[-1]],
            "mean_information_content": _r(
                statistics.mean(sites[i - 1]["information_content"] for i in span), 3
            ),
            "mean_gap_fraction": _r(statistics.mean(sites[i - 1]["gap_fraction"] for i in span), 3),
        }
    occupied = [x for x in sites if x["gap_fraction"] <= 0.5]
    occ_ic = [x["information_content"] for x in occupied if x["site_index"] not in motif_sites]
    rank = sorted(occupied, key=lambda x: -(x["information_content"] or 0))
    top20 = [x["site_index"] for x in rank[:20]]

    # Physicochemical variability: SD of Kyte-Doolittle hydropathy among residues per column.
    pc = physicochemical(profiler.alignment, "protein").array
    resid = (pc[:, :, 5] == 0) & (pc[:, :, 6] == 0)
    hyd = np.where(resid, pc[:, :, 0] * 4.5, np.nan)
    sd = np.nanstd(hyd, axis=0)
    occ_idx = np.array([x["site_index"] - 1 for x in occupied])
    motif_idx = np.array(sorted(i - 1 for i in motif_sites))
    other_idx = np.setdiff1d(occ_idx, motif_idx)

    # Figure panels.
    PANELS.mkdir(parents=True, exist_ok=True)
    profiler.plot_alignment(
        PANELS / "ras_alignment.png",
        start=1,
        end=50,
        max_rows=20,
        order="consensus_identity",
        tracks=("conservation", "gaps"),
        dpi=300,
    )
    profiler.plot_logo(PANELS / "ras_logo_g1.png", start=1, end=20, dpi=300)
    profiler.to_image(PANELS / "ras_rgb.png", kind="rgb", scale=4)
    profiler.to_image(PANELS / "ras_labels.png", kind="label")

    # Encodings: exact round trip and a windowed export consumed by PyTorch.
    cat = categorical(profiler.alignment, "protein")
    decoded = decode_categorical(cat.array, cat.ids, full_vocabulary("protein"))
    roundtrip = [r.sequence for r in decoded.records] == [
        r.sequence.upper() for r in profiler.alignment.records
    ]
    t1 = time.perf_counter()
    profiler.export_encoding(
        RESULTS / "ras_ml", encoding="onehot", window=64, stride=32, rows_per_tile=32, stem="ras"
    )
    export_seconds = time.perf_counter() - t1
    npz = next((RESULTS / "ras_ml").glob("ras_onehot.npz"))
    ds = AlignmentTensorDataset(npz, channels_first=True)
    import torch
    from torch.utils.data import DataLoader

    loader = DataLoader(
        ds,
        batch_size=8,
        shuffle=False,
        collate_fn=lambda b: {k: torch.stack([to_torch(i[k]) for i in b]) for k in b[0]},
    )
    first = next(iter(loader))
    manifest = json.loads(npz.with_name("ras_onehot.manifest.json").read_text())
    return {
        "dataset": "Pfam PF00071 (Ras family) seed alignment",
        "stats_seconds_exact_pairwise": round(stats_seconds, 3),
        "sequences": s["number_of_sequences"],
        "sites": s["alignment_length"],
        "gap_fraction": _r(s["gap_fraction"]),
        "high_gap_sites_gt_50pct": sum(1 for x in sites if x["gap_fraction"] > 0.5),
        "mean_pairwise_identity": _r(s["mean_pairwise_identity"]),
        "pairwise_identity_range": [
            _r(s.get("min_pairwise_identity")),
            _r(s.get("max_pairwise_identity")),
        ],
        "qc_findings": [
            {"code": f.code, "severity": f.severity, "count": f.count}
            for f in profiler.qc_findings()
        ],
        "g_box_motifs": motifs,
        "occupied_sites": len(occupied),
        "median_information_content_other_occupied": _r(statistics.median(occ_ic), 3),
        "motif_sites_in_top20_information": sum(1 for i in top20 if i in motif_sites),
        "motif_sites_total": len(motif_sites),
        "hydropathy_sd_median_motif": _r(float(np.nanmedian(sd[motif_idx])), 3),
        "hydropathy_sd_median_other_occupied": _r(float(np.nanmedian(sd[other_idx])), 3),
        "categorical_roundtrip_exact": bool(roundtrip),
        "onehot_windows": {
            "tiles": len(ds),
            "tile_shape_channels_first": list(first["x"].shape[1:]),
            "first_batch_shape": list(first["x"].shape),
            "mask_shape": list(first["mask"].shape),
            "export_seconds": round(export_seconds, 3),
            "channels": manifest.get("channels"),
            "manifest_keys": sorted(manifest),
        },
    }


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    out = {"influenza_h3n2_ha": influenza(), "ras_pf00071": ras()}
    (RESULTS / "use_cases.json").write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
