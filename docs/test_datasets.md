# Test Datasets

MSATK includes small synthetic test alignments under `tests/data/`.

These datasets are designed to test format parsing, molecule-type detection, summary statistics, codon-aware analysis, protein summaries, and edge-case behavior.

## DNA Datasets

- `tiny_dna.fasta`: small aligned DNA example
- `dna_with_gaps.fasta`: tests gap counting and site missingness
- `dna_with_ambiguous_bases.fasta`: tests ambiguous character handling
- `unequal_lengths.fasta`: tests strict validation failure

## Codon Datasets

- `valid_codon_alignment.fasta`: valid CDS alignment
- `codon_with_stop.fasta`: contains an internal stop codon
- `codon_not_divisible_by_three.fasta`: invalid codon-mode input
- `codon_with_frameshift.fasta`: contains a gap inside a codon-like alignment

## Protein Datasets

- `tiny_protein.faa`: small protein alignment
- `protein_with_gaps.faa`: tests protein gap handling
- `protein_with_unknowns.faa`: tests unknown residue handling

## Format Datasets

- `tiny.aligned.fasta`
- `tiny.a3m`
- `tiny.phy`
- `tiny_relaxed.phylip`
- `tiny.clustal`
- `tiny.nexus`
- `tiny.stockholm`
- `tiny.maf`
- `tiny.sam`
- `tiny.bam`
- `tiny.cram`: generated when `pysam`/HTSlib is available
- `reference.fasta`: reference used to generate CRAM

The BAM and CRAM fixtures are generated from the synthetic reads in `tiny.sam` using `create_binary_alignment_fixtures.py`. They require `pysam` for regeneration and parsing.

## Edge Cases

- `empty_file.fasta`
- `single_sequence.fasta`
- `all_gap_column.fasta`
- `duplicate_ids.fasta`
- `invalid_characters.fasta`
- `very_long_sequence_names.fasta`
