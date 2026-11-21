# MSATK logo

The current logo is **v4** (2026-11-21). It is drawn by `make_logo.py`, so it can be
regenerated or changed exactly:

```bash
python logo/make_logo.py
```

| File | Use |
| --- | --- |
| `msatk_logo.svg` / `.png` | Full logo for light backgrounds (README, docs pages) |
| `msatk_logo_dark.svg` / `.png` | Full logo for dark backgrounds |
| `msatk_mark.svg` / `.png` | Square mark: docs header, favicon, HTML report |
| `msatk_mark_dark.svg` / `.png` | Square mark for dark backgrounds |

All files have transparent backgrounds. SVG text is converted to paths, so no fonts are
needed to display them. The PNG copies used elsewhere are
`docs/assets/msatk_logo.png`, `docs/assets/msatk_mark.png` and
`src/msatk/assets/msatk_logo.png` (the mark, used by `msatk profile` reports). After
regenerating, copy them again.

The mark is a real five-sequence alignment next to its guide tree. Cells are coloured by
consensus support: teal where every sequence agrees, light teal where most agree, and
grey for substitutions and weakly supported columns. Dashes are gaps.

## Previous versions

Kept in `archive/` for reference; do not use them in new material.

| Version | File | Notes |
| --- | --- | --- |
| v1 | `archive/msatk_logo_v1.png` | Navy tree and alignment with wordmark; used in docs and reports until v0.2.0 |
| v1 (green) | `archive/msatk_logo_v1_green.png` | Green colourway of v1 |
| v2 | `archive/msatk_logo_v2.png` | Blue DNA helix (raster, fake transparency baked in) |
| v3 | `archive/msatk_logo_v3.png` | Blue helix with glow on grey; README logo until v0.2.0 |
