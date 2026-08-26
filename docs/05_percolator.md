# 5. Filtering matches with Percolator

Percolator is a semi-supervised machine-learning tool that separates real
peptide identifications from noise, using the target/decoy matches from
tide-search to learn what "real" vs. "random" scores look like, then
reports only matches below a false discovery rate (FDR) threshold you
choose.

You do **not** need to install Percolator separately — it's bundled inside
crux and invoked as `crux percolator`.

## 5.1 Run it

```bash
cd proteomics_pipeline_course/mac_linux/scripts   # or windows/scripts on Windows
./04_percolator.sh /path/to/search_out
```

This reads `search_out/tide-search.txt` and writes several files into the
same folder, most importantly:

- **`percolator.target.psms.txt`** — every target PSM, each with a
  computed `q-value` column. **Despite `--test-fdr`/`--train-fdr` being
  set, this file is NOT pre-filtered down to just the confident matches**
  — percolator uses those FDR values to *train* its scoring model, but
  still writes out every target PSM regardless of where its q-value ends
  up. **This is the file the next step (moFF) needs** — but the actual
  FDR filtering happens in that next step's script
  (`05_percolator_to_moff.py`'s `--max-q-value`, see doc 6), not here.
- `percolator.target.peptides.txt` — same idea, collapsed to one row per
  unique peptide rather than per PSM.
- `percolator.decoy.*.txt` — the decoy-side matches (useful for diagnostics,
  not for biological interpretation).

Typically finishes in well under a minute for a modest dataset, once
tide-search has finished.

## 5.2 Understanding the output columns

| Column | Meaning |
|--------|---------|
| `PSMId` | internal ID, encodes which input file / scan / charge / rank this came from |
| `filename` | which mzML file this match came from |
| `score` | percolator's learned discriminant score (higher = more confident) |
| `q-value` | the FDR threshold at which this specific match would first be kept — e.g. `q-value = 0.03` means this match is part of the set kept at 3% FDR or looser |
| `posterior_error_prob` | estimated probability *this individual match* is wrong |
| `peptide` | the identified peptide sequence (flanking residues + any modification notation) |
| `proteinIds` | protein(s) this peptide maps to |
| `rt` | retention time (minutes) — the converter script in doc 6 does not actually rely on this column (it re-derives retention time from the mzML file directly), but it's there if you want it |

## 5.3 Choosing an FDR threshold

`04_percolator.sh` uses `--test-fdr 0.1 --train-fdr 0.1` (10%) to *train*
percolator's model, matching this lab's established defaults for
exploratory/QC-level searches on small databases — this mostly affects how
percolator learns to separate real from decoy matches, not what ends up in
the output file (see doc 5.1: that file isn't pre-filtered regardless of
this setting). **The threshold that actually determines what counts as
"confident" is `--max-q-value` on `05_percolator_to_moff.py` in the next
step** — it defaults to matching `0.1` here, but **for a real study you
will usually want a stricter threshold — 1% FDR (`--max-q-value 0.01`) is
the standard reported in publications.** Ask your instructor/PI which
threshold is expected for your project. You can change `--max-q-value`
freely without rerunning percolator (it's just filtering the existing
q-value column) — but if you also want percolator's own model retrained at
a different FDR, that does require rerunning `04_percolator.sh` with
matching `--test-fdr`/`--train-fdr`.

## 5.4 A quick quality check

Before moving on, it's worth checking the search found something
reasonable. Note `percolator.target.psms.txt` itself is **not** filtered
(see doc 5.1) — for a real confident-PSM count, filter by q-value first:

```bash
wc -l search_out/percolator.target.psms.txt                          # total target PSMs (NOT all confident)
awk -F'\t' 'NR>1 && $4<=0.1' search_out/percolator.target.psms.txt | wc -l  # confident at 10% FDR (column 4 = q-value)
cut -f7 search_out/percolator.target.psms.txt | sort -u | wc -l       # roughly how many distinct proteins (unfiltered)
```

For a real sample against an appropriately sized database, expect at least
tens to hundreds of confident PSMs and more than a handful of distinct
proteins. If you get suspiciously few (e.g. single digits), double check:
the right database was used, the FASTA matches the organism, and the mzML
files converted correctly (doc 3).

Next: [6. Quantification with moFF](06_spectral_count_moff.md)
