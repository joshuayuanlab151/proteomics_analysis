# 6. Quantification with moFF

Percolator (doc 5) told you *which* peptides were confidently identified.
moFF answers the next question: *how much* of each peptide was present in
each sample. It does this by going back to the original mzML signal and
measuring the actual MS1 peak intensity for every confidently identified
peptide, in every run — including, optionally, in runs where that specific
peptide wasn't identified by MS2 but likely is still present (see 6.3).

## 6.1 Step A: convert percolator output into moFF's input format

moFF needs a specific tab-delimited format (`peptide`, `prot`, `mod_peptide`,
`rt`, `mz`, `mass`, `charge`) — crux's percolator output isn't in that
format, and (depending on your crux version) may or may not include a
usable retention time column, but never includes m/z or mass. A converter
script bridges the two automatically:

```bash
conda activate proteomics    # needs lxml, which lives in this env, not moff
python 05_percolator_to_moff.py \
    --percolator search_out/percolator.target.psms.txt \
    --tide-search search_out/tide-search.txt \
    --mzml-dir /path/to/mzML \
    --out-dir moff_input
```

(Run this from inside `mac_linux/scripts/` or `windows/scripts/`.)

This writes one file per run into `moff_input/` (e.g. `1.txt`, `2.txt`),
named to match each run's mzML filename — moFF requires that exact naming
convention. Retention times are read directly out of each mzML file by
scan number rather than trusted from crux's output, so this works
regardless of crux version.

## 6.2 Step B: run moFF

```bash
conda activate moff           # different environment! moFF needs old pinned packages
./06_run_moff.sh moff_input /path/to/mzML moff_out project
```

For each peptide in each run's input file, moFF extracts the extracted-ion
chromatogram (XIC) around that peptide's identified retention time and m/z,
finds the peak apex, and reports its intensity. Once installed correctly
(indexed mzML + the bugfix patch described below), the full pipeline on
the bundled 2-run example dataset — including RAW conversion, search,
percolator, and moFF quantification with match-between-runs on — normally
completes in a few minutes on an otherwise-idle machine. It can take much
longer if the machine is busy with other work at the same time (see doc 7)
— that's a system-load effect, not a sign of a problem. For a real, larger
study with many samples, expect this step to take longer — let it run in
the background if needed.

**Note:** `00_setup_environment.sh` applies a small patch to moFF
(`patches/moff_random_access_fix.patch`) that fixes a real bug in the
unmaintained upstream tool — without it, moFF silently reports zero peaks
for every single peptide instead of erroring. See
[doc 7](07_troubleshooting.md) if you ever see that happen (e.g. after a
manual, non-scripted install).

## 6.3 What comes out

- One output file per input run (e.g. `1_moff_result.txt`), with all the
  original columns plus `intensity`, `rt_peak` (the actual apex RT found),
  `SNR`, and a `matched` flag (`0` = confirmed by MS2 identification in
  this run, `1` = quantified only via match-between-runs, see below).
- **`peptide_summary_intensity_project.tab`** — the file you'll actually
  use for downstream analysis: one row per peptide, one column per run,
  values are intensity. This is your quantification table, ready to import
  into R/Python/Excel for statistics (t-tests, fold-change, PCA, etc.).

### Match-between-runs (MBR)

`06_run_moff.sh` runs with `--mbr on` by default, matching this lab's
established workflow. This means: if peptide X was confidently identified
by MS2 in run A but missed in run B (e.g. just below the detection
threshold that scan), moFF will still look for and report its MS1
intensity in run B, using the retention-time relationship learned between
runs A and B. This recovers more usable data points, at the cost of some
values in your table coming from inference rather than a direct
identification in that specific run — check the `matched` column if you
need to know which is which. Set `--mbr off` in the script (or pass a
different value when calling `moff_all.py` directly) if your analysis
requires only directly-identified values.

## 6.4 Sanity-check the result

Before treating the numbers as final, check that the run actually worked:

```bash
wc -l moff_out/peptide_summary_intensity_project.tab
head moff_out/peptide_summary_intensity_project.tab
```

You should see one row per peptide and reasonable-looking (non-zero,
non-identical-across-all-rows) intensity values. If every intensity is 0 or
missing, the most common cause is a retention-time or m/z mismatch between
the identification file and the mzML — double check that step A used the
mzML files from the *same* conversion run as step 3/4/5 (not a stale or
mismatched copy).

Next: [7. Troubleshooting](07_troubleshooting.md)
