# Application Note analyses

Scripts that regenerate every number, table and figure in `paper/paper.md` and
`paper/supplement.md`.

```bash
python paper/analysis/prepare_data.py            # download Pfam seeds + H3N2 HA CDS (or --offline)
python paper/analysis/compare_tools.py --amas-python /path/to/python-with-amas
python paper/analysis/use_cases.py               # both use cases + figure panels
python paper/analysis/make_figure1.py            # paper/figures/figure1.{png,pdf}
```

`compare_tools.py` needs AMAS (`pip install amas`) in the interpreter you pass;
it is kept separate so that MSATK's own environment is unchanged. `use_cases.py`
needs PyTorch for the `DataLoader` check.

| Path | Contents |
|---|---|
| `data/` | Input alignments (`*.fasta`), raw downloads (`*.sto.gz`) and H3N2 host labels (`h3n2_ha_metadata.tsv`) |
| `results/tool_comparison.{csv,json}` | Table 2 / Table S1 |
| `results/use_cases.json` | Every number quoted in the use cases (Tables S3–S6) |
| `results/panels/` | MSATK-rendered panels used in Figure 1 |
| `results/h3n2/`, `results/ras_ml/` | Full MSATK outputs (regenerated, not committed) |
