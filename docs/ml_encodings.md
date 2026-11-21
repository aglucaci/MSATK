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
| `rgb` | `(N, L, 3)` uint8 | Theme colours for RGB-pretrained models; **lossy** (colours are labels) |

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

## Alignments of different sizes: export strategies

Alignments differ in size (40 × 300, 120 × 900, 300 × 1500 sequences × sites), while
many models batch fixed-size inputs. MSATK keeps the alignment intact and offers three
strategies. None of them resamples: every exported cell is one alignment cell.

| Strategy | Output `x` | Use it for | Information kept |
|---|---|---|---|
| `native` | `(N, L, C)` | Models that accept variable sizes (CNNs with adaptive pooling, custom architectures) | Everything |
| `tiles` | `(T, height, width, C)` | CNNs and ViTs with a fixed input size; long alignments (the recommended default for ML) | Everything: tiles cover every cell and record their coordinates |
| `pad` | `(1, height, width, C)` | Collections of alignments that all fit one target size | Everything, if the alignment fits; larger inputs are rejected unless `--overflow crop` (lossy, recorded) |

Resizing to a fixed size such as 224 × 224 merges or duplicates residues and is not
offered as a strategy; `msatk image --size` exists for display only and is marked lossy.

```bash
# Recommended: tiles that keep every cell (last tile on each axis padded)
msatk encode alignment.fasta --strategy tiles --height 128 --width 256 --out ml/

# Overlapping tiles
msatk encode alignment.fasta --strategy tiles --width 256 --stride 128 --out ml/

# One fixed-size tensor; fails if the alignment is larger than the target
msatk encode alignment.fasta --strategy pad --height 512 --width 2048 --out ml/
```

```python
msa.export_encoding("ml/", encoding="onehot", strategy="tiles", height=128, width=256,
                    labels={"split": "train"})
```

`tiles` and `pad` NPZ files hold `x`, `mask` (True = real cell, False = padding) and
`tiles` (`row_start, row_end` 0-based in the exported row order; `site_start, site_end`
1-based in the original alignment), so the full alignment can be reassembled exactly.
The manifest's `strategy` block records the target size, strides, the number of padded
cells and any cropping. (`window=`/`--window` and `rows_per_tile=` still work as
older names for `width` and `height`.)

### Codon alignments

For `--type codon`, tile `width` and `stride` must be multiples of 3 and `--start` must
be codon position 1, so every tile starts on a codon boundary. MSATK refuses settings
that would break the reading frame, and records `reading_frame` in the manifest.

### Row order

The order of sequences changes the spatial patterns a CNN sees, so it is recorded with
the export (`row_order` in the manifest; ordering is computed on the full alignment):

| `--order` | Order |
|---|---|
| `input` (default) | File order |
| `reference` | Decreasing identity to `--order-reference ID` (default: first sequence) |
| `consensus_identity` | Decreasing identity to the consensus |
| `gaps` | Increasing gap fraction |
| `cluster` | Leaf order of average-linkage (UPGMA) clustering on p-distance; deterministic, pure NumPy |
| `file:PATH` | IDs listed one per line; unlisted IDs follow in file order |

### Pretrained RGB models

Pretrained vision models expect three channels, whereas one-hot DNA has six or more.
Two options:

* `--encoding rgb` exports `(…, 3)` uint8 tiles coloured by a theme. Colours are display
  labels, so this is marked lossy; use it for compatibility, not as the scientific
  representation.
* Keep the one-hot channels and replace the model's first convolution with one that
  takes `C` input channels.

### Datasets and batching

```python
from msatk.ml import AlignmentTensorDataset, pad_batch, to_torch

# Several alignments exported with the same --height/--width form one dataset.
ds = AlignmentTensorDataset(["ml/a_onehot.npz", "ml/b_onehot.npz"], channels_first=True)
item = ds[0]                    # {"x", "mask", "tile", "source"}

# Native exports of different sizes: pad each batch to its largest member.
native = AlignmentTensorDataset(["n/a_onehot.npz", "n/b_onehot.npz"], channels_first=True)
batch = pad_batch([native[0], native[1]], multiple_of=16)   # x, mask, source
```

`AlignmentTensorDataset` implements `__len__`/`__getitem__`, so it works directly with
`torch.utils.data.DataLoader` (use `collate_fn=pad_batch` for native exports).
`to_torch` and `to_tensorflow` convert arrays.

## Images

```bash
msatk image alignment.fasta --kind rgb --scale 4 --out aln.png      # lossless integer upscale
msatk image alignment.fasta --kind label --out aln_labels.png       # pixel = vocabulary index
msatk image alignment.fasta --size 224x224 --out aln_224.png        # lossy, recorded in manifest
```

Resizing (for RGB-pretrained models) is nearest-neighbour and is recorded as a
lossy transform; prefer `msatk encode --strategy tiles` for scientific work.

## Transforms

`msatk.encoding.subset_alignment` crops sites, selects/reorders rows and drops
gap-rich columns while returning a `CoordinateMap` from output to original rows
and sites.
