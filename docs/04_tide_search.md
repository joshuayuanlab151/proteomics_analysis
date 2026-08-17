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
together is fine and normal (with `--concat T`, they all end up in a
single combined `tide-search.txt`, tagged by which file each match came
from).

**This is the slowest step of the pipeline.** For a database with a few
thousand proteins and several runs, budget anywhere from tens of minutes
to a few hours depending on your computer and how many spectra you're
searching. Let it run in the background if needed — there's no need to
babysit it.

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
