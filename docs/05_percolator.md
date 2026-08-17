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

- **`percolator.target.psms.txt`** — the filtered, confident peptide
  identifications. **This is the file the next step (moFF) needs.**
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

`04_percolator.sh` uses `--test-fdr 0.1 --train-fdr 0.1` (10%), matching
this lab's established defaults for exploratory/QC-level searches on small
databases. **For a real study you will usually want a stricter threshold —
1% FDR (`--test-fdr 0.01 --train-fdr 0.01`) is the standard reported in
publications.** Ask your instructor/PI which threshold is expected for your
project, and edit the script (or copy it) accordingly. Filtering by
`q-value` at a stricter threshold, without rerunning percolator, is also
possible by filtering the output table with e.g. `awk` or in R/Python — but
rerunning with a matching `--test-fdr`/`--train-fdr` is the recommended way.

## 5.4 A quick quality check

Before moving on, it's worth confirming the search found something
reasonable:

```bash
wc -l search_out/percolator.target.psms.txt      # how many confident PSMs?
cut -f7 search_out/percolator.target.psms.txt | sort -u | wc -l   # roughly how many distinct proteins?
```

For a real sample against an appropriately sized database, expect at least
tens to hundreds of confident PSMs and more than a handful of distinct
proteins. If you get suspiciously few (e.g. single digits), double check:
the right database was used, the FASTA matches the organism, and the mzML
files converted correctly (doc 3).

Next: [6. Quantification with moFF](06_spectral_count_moff.md)
