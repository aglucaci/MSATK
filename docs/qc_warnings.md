# QC findings

MSATK evaluates explicit rules and reports **severity-ranked findings** (error,
warning, info, ok), each with a code, scope, count, explanation, suggested follow-up
and the affected sequences or sites. Findings appear in the terminal,
`qc_findings.{json,csv}`, `qc_warnings.txt` and the report.

| Code | Severity | Trigger |
|---|---|---|
| `UNEQUAL_LENGTHS` | warning | Input rows differ in length (padded as missing data) |
| `DUPLICATE_IDS` | warning | Identifier occurs more than once |
| `DUPLICATE_SEQUENCES` | info | Identical sequences |
| `INVALID_CHARACTERS` | warning | Characters outside the molecule alphabet |
| `HIGH_GAP_SEQUENCES` | warning | Sequence gap fraction > `--max-gap-seq` (0.3) |
| `HIGH_GAP_SITES` | warning | Column gap fraction > `--max-gap-site` (0.5) |
| `UNOBSERVED_SITES` | info | Columns with no canonical residue |
| `DIVERGENT_SEQUENCES` | warning | Consensus identity < 0.5 or robust z < −3.5 |
| `COMPOSITION_OUTLIERS` | warning | Composition-distance robust z > 3.5 |
| `HIGH_AMBIGUITY_SEQUENCES` | warning | > 5 % ambiguity codes |
| `NON_OVERLAPPING_PAIRS` | warning | Pairs sharing no comparable site |
| `CODON_LENGTH_NOT_MULTIPLE_OF_3` | warning | Codon mode, L mod 3 ≠ 0 |
| `INTERNAL_STOP_CODON` | warning | Codon mode, in-frame stop before the last codon |
| `PARTIAL_CODON_GAPS` | warning | Codon mode, gaps that split codons |
| `INTERNAL_STOPS_PROTEIN` | warning | Protein `*` before the sequence end |
| `SINGLE_SEQUENCE` | info | Only one sequence |
| `NO_MAJOR_ISSUES` | ok | No warning triggered |

Thresholds are adjustable from the CLI and via `msatk.core.qc.QCThresholds`.
