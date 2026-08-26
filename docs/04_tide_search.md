# 4. Peptide identification with crux tide-search

## 4.1 Build the search index (once per database)

```bash
cd proteomics_pipeline_course/mac_linux/scripts   # or windows/scripts on Windows
./02_build_database_index.sh /path/to/database.fasta /path/to/tide_index
```

(No `conda activate` needed for this step — the script finds crux directly
under `software/crux/`, installed by doc 1 / 1b's setup script.)

What this does (see the script for the exact command): it digests every
protein into peptides (using trypsin cleavage rules by default), generates
a matching decoy peptide for each target peptide by reversing it, and
writes an indexed database that tide-search can search against quickly.
`--mods-spec C+57.02146` tells it every cysteine (C) carries a fixed
+57.02146 Da modification (carbamidomethylation) — the near-universal
result of the iodoacetamide alkylation step in standard sample prep. If
your samples were prepared differently, ask your instructor whether this
should change.

You only need to rebuild the index if your FASTA changes — reuse the same
index for every sample searched against that database.

## 4.2 Run the search

```bash
./03_tide_search.sh /path/to/tide_index /path/to/search_out /path/to/mzML/*.mzML
```

You can pass one or many `.mzML` files at once — passing several samples
together is fine and normal. There's no special flag for this: crux's
own usage text describes it as `crux tide-search [options] <tide spectra
file> <tide database>`, where *"multiple files can be included on the
command line (space delimited), prior to the name of the database"* —
that's exactly what `/path/to/mzML/*.mzML` above does (the shell expands
the glob into however many space-separated filenames match, before
`03_tide_search.sh` appends the database path after them). With
`--concat T`, they all end up in a single combined `tide-search.txt`,
tagged by which file each match came from — `03_tide_search.sh` sets
this by default.

**This is the slowest step of the pipeline.** For a database with a few
thousand proteins and several runs, budget anywhere from tens of minutes
to a few hours depending on your computer and how many spectra you're
searching. Let it run in the background if needed — there's no need to
babysit it.

## 4.2.1 `--mz-bin-width`: match this to your instrument's resolution

`03_tide_search.sh` currently hardcodes `--mz-bin-width 1.0005079`.
**This is only correct for low-resolution ion trap MS2 data (CID
fragmentation).** Straight from crux's own `--help` text: *"For low
resolution ion trap ms/ms data 1.0005079 and for high resolution ms/ms
0.02 is recommended. Default = 0.02."* — note that crux's actual default
is the **high-res** value, not `1.0005079`; easy to get backwards (this
doc used to). `1.0005079` is used here because this pipeline's own
bundled example data really is low-resolution ion trap — verified, not
assumed. **If your own data is high-resolution (Orbitrap/HCD fragmentation),
edit the script to use `--mz-bin-width 0.02` instead** — using the
low-res bin width on high-res data throws away real discriminating power
in the search.

**Check which one your data actually is, rather than guessing:**
ThermoRawFileParser embeds the instrument's own per-scan description
(cvParam `MS:1000512`, "filter string") in the mzML. `FTMS` = Orbitrap/FT
(high-res); `ITMS` = ion trap (low-res):
```bash
grep -o 'accession="MS:1000512" value="[^"]*"' yourfile.mzML \
    | sed 's/.*value="//;s/"$//' | cut -c1-4 | sort | uniq -c
```
Some methods deliberately mix modes (high-res MS1 survey + fast low-res
MS2) — check what your MS2 scans specifically are, since that's what
`--mz-bin-width` actually discretizes. Crux can also estimate this itself
instead of you setting it by hand: `--auto-mz-bin-width warn` (falls back
to the default if estimation fails) or `fail` (errors instead of silently
guessing wrong).

## 4.3 What comes out

`search_out/tide-search.txt` — one row per candidate peptide-spectrum
match. Important columns:

| Column | Meaning |
|--------|---------|
| `scan` | which spectrum in the mzML this match is for |
| `charge` | the peptide's inferred charge state |
| `retention time` | when this spectrum was measured during the run (minutes) |
| `spectrum precursor m/z` | measured mass-to-charge of the precursor ion |
| `peptide mass` | the candidate peptide's theoretical mass |
| `xcorr score` | the main match-quality score (higher = better) |
| `sequence` | the candidate peptide sequence (with flanking residues, e.g. `K.PEPTIDE.R`) |
| `protein id` | which protein(s) this peptide could come from |
| `target/decoy` | whether this match is against a real ("target") or fake ("decoy") sequence |

**This file is NOT your final answer.** At this stage the file contains
many matches that pass simply by chance — roughly comparable numbers of
target and decoy top-scoring matches would appear if none of your peptides
were real, which is exactly why the next step (percolator) exists.

Next: [5. Filtering with percolator](05_percolator.md)
