"""Profile hidden Markov models: HMMER3 (``.hmm``) and HH-suite (``.hhm``).

MSATK reads model metadata and, for HMMER3, match-state emission probabilities,
from which it reports per-state emission entropy and information content relative
to the model's own mean composition (``COMPO``). Profile-to-sequence scoring is
not performed.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from statistics import mean
from typing import Any

from msatk.exceptions import AlignmentFormatError


@dataclass
class ProfileModel:
    kind = "profile_model"
    fmt: str
    header: dict[str, str]
    alphabet: list[str] = field(default_factory=list)
    match_emissions: list[list[float]] = field(default_factory=list)  # probabilities
    composition: list[float] | None = None
    source: str = ""

    def state_rows(self) -> list[dict[str, Any]]:
        rows = []
        for k, probs in enumerate(self.match_emissions, start=1):
            ent = -sum(p * math.log2(p) for p in probs if p > 0)
            ic = None
            if self.composition:
                ic = sum(
                    p * math.log2(p / q) for p, q in zip(probs, self.composition) if p > 0 and q > 0
                )
            best = max(range(len(probs)), key=probs.__getitem__) if probs else None
            rows.append(
                {
                    "match_state": k,
                    "consensus": self.alphabet[best] if best is not None else None,
                    "max_emission": max(probs) if probs else None,
                    "entropy_bits": ent,
                    "relative_entropy_bits": ic,
                }
            )
        return rows

    def summary(self) -> dict[str, Any]:
        h = self.header
        rows = self.state_rows()

        def num(key: str) -> float | None:
            match = re.match(r"\s*([-+]?\d+(?:\.\d*)?)", h.get(key, ""))
            return float(match.group(1)) if match else None

        return {
            "alignment_kind": "profile_model",
            "detected_format": self.fmt,
            "model_name": h.get("NAME"),
            "model_accession": h.get("ACC"),
            "alphabet": h.get("ALPH"),
            "model_length": int(num("LENG") or 0) or None,
            "training_sequences": int(num("NSEQ") or 0) or None,
            "effective_sequences": num("EFFN") or num("NEFF"),
            "mean_match_entropy_bits": mean(r["entropy_bits"] for r in rows) if rows else None,
            "mean_match_relative_entropy_bits": mean(
                r["relative_entropy_bits"] for r in rows if r["relative_entropy_bits"] is not None
            )
            if rows and self.composition
            else None,
        }

    def tables(self) -> dict[str, list[dict[str, Any]]]:
        return {"match_states": self.state_rows()}


def _prob(token: str) -> float:
    return 0.0 if token == "*" else math.exp(-float(token))


def parse_hmmer(text: str, source: str = "") -> ProfileModel:
    lines = text.splitlines()
    if not lines or not lines[0].startswith("HMMER3"):
        raise AlignmentFormatError("Not a HMMER3 profile (expected a 'HMMER3' first line).")
    header: dict[str, str] = {}
    model = ProfileModel("hmmer", header, source=source)
    i = 1
    while i < len(lines) and not lines[i].startswith("HMM "):
        parts = lines[i].split(None, 1)
        if len(parts) == 2:
            header.setdefault(parts[0], parts[1].strip())
        i += 1
    if i >= len(lines):
        raise AlignmentFormatError("HMMER3 file has no 'HMM' section.")
    model.alphabet = lines[i].split()[1:]
    k = len(model.alphabet)
    i += 2  # alphabet line, transition header line
    for line in lines[i:]:
        parts = line.split()
        if not parts or parts[0] == "//":
            break
        if parts[0] == "COMPO" and len(parts) >= k + 1:
            model.composition = [_prob(t) for t in parts[1 : k + 1]]
        elif parts[0].isdigit() and len(parts) >= k + 1:
            model.match_emissions.append([_prob(t) for t in parts[1 : k + 1]])
    return model


def parse_hhm(text: str, source: str = "") -> ProfileModel:
    header: dict[str, str] = {}
    for line in text.splitlines():
        if line.startswith("#") or line.startswith("HMM"):
            break
        parts = line.split(None, 1)
        if len(parts) == 2:
            header.setdefault(parts[0], parts[1].strip())
    if "LENG" not in header:
        raise AlignmentFormatError("Not an HH-suite HHM profile (no LENG line).")
    return ProfileModel("hhm", header, source=source)
