# Machine-learning representations

MSATK separates **human visualisation** (RGB colours are arbitrary labels) from
**faithful numerical encodings**. Every export writes a `manifest.json` describing
row order, sequence IDs, original site coordinates, vocabulary, channel semantics,
dtype, padding, transforms, labels, input SHA-256 and MSATK version.

## Encodings

| Encoding | Shape | Notes |
|---|---|---|
| `categorical` | `(N, L)` uint8 | Full vocabulary (IUPAC / extended amino acids, `-`, `.`, `?`, `*`); **exactly reversible** |
| `onehot` | `(N, L, K+2)` or `(N, L, K+3)` float32 | Canonical residues + `<gap>` + `<missing>` (+ `<ambiguous>`) |
| `physicochemical` | `(N, L, 7)` float32 | Proteins: hydropathy, charge, volume, polar, aromatic, gap mask, non-residue mask |
| `features` | `(N, L, 7)` float32 | residue/gap/missing/ambiguous masks, matches-consensus, site normalised entropy, site occupancy |

One-hot ambiguity policies (`--ambiguity`):

* `missing` — ambiguity codes go to the `<missing>` channel;
* `fractional` — IUPAC codes are spread over compatible bases (`R` → A 0.5, G 0.5;
  `N` → 0.25 each); `?` stays missing;
* `separate` — ambiguity codes get an `<ambiguous>` channel.

Every non-padding cell sums to exactly 1 across one-hot channels.

```python
from msatk import MSATK
from msatk.encoding import decode_categorical, full_vocabulary

msa = MSATK("alignment.fasta")
t = msa.to_tensor("categorical")
assert decode_categorical(t.array, t.ids, full_vocabulary("dna")).sequences == \
    [s.upper() for s in msa.alignment.sequences]

x = msa.to_tensor("onehot", ambiguity="fractional").array   # numpy (N, L, 6)
```

## Windows and tiles

Biological alignments are not photographs: rotation, mirroring, interpolation and
colour jitter destroy their meaning. MSATK provides integer-exact windowing with
padding masks instead of resizing:

```python
msa.export_encoding("ml/", encoding="onehot", window=256, stride=128, rows_per_tile=64,
                    labels={"split": "train"})
```

The NPZ holds `x` `(tiles, rows, width, channels)`, `mask` (True = real cell) and
`tiles` (`row_start, row_end, site_start, site_end` in original coordinates).

```python
from msatk.ml import AlignmentTensorDataset, to_torch

ds = AlignmentTensorDataset("ml/alignment_onehot.npz", channels_first=True)
item = ds[0]                    # {"x": ..., "mask": ..., "tile": ...}
batch = to_torch(item["x"])     # optional PyTorch adapter
```

`AlignmentTensorDataset` implements `__len__`/`__getitem__`, so it works directly
with `torch.utils.data.DataLoader`. `to_tensorflow` is the TensorFlow adapter.

## Images

```bash
msatk image alignment.fasta --kind rgb --scale 4 --out aln.png      # lossless integer upscale
msatk image alignment.fasta --kind label --out aln_labels.png       # pixel = vocabulary index
msatk image alignment.fasta --size 224x224 --out aln_224.png        # lossy, recorded in manifest
```

Resizing (for RGB-pretrained models) is nearest-neighbour and is recorded as a
lossy transform; prefer windows for scientific work.

## Transforms

`msatk.encoding.subset_alignment` crops sites, selects/reorders rows and drops
gap-rich columns while returning a `CoordinateMap` from output to original rows
and sites.
