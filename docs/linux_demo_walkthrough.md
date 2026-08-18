# Linux command-line walkthrough: RAW → final protein/gene list

A back-to-back command reference for demoing the identification pipeline
by hand on Linux: CD-HIT → build the search database → ThermoRawFileParser
→ crux tide-search → percolator → protein/gene list → spectral counts.
Every command and log excerpt below was actually run against real test
data (`PSM_BL.fasta`, one real `.raw` file) while writing this doc — not
guessed. MS1-intensity-based quantification (moFF) is a separate, later
topic, not covered here.

## Prerequisites

Assumes crux, ThermoRawFileParser, and a conda environment with `lxml`
are already installed, and that you're working inside that conda
environment. CD-HIT is not part of that environment by default — install
it once:

```bash
conda activate proteomics
conda install -y -c bioconda -c conda-forge cd-hit
```

Set these once per session so every command below can just be copy-pasted:

```bash
export CRUX=/path/to/mac_linux/software/crux/bin/crux
export TRFP=/path/to/mac_linux/software/ThermoRawFileParser/ThermoRawFileParser
export CONTAM=/path/to/reference_data/common_contaminants.fasta
export WORKDIR=~/demo_run
mkdir -p "$WORKDIR" && cd "$WORKDIR"
```

---

## Step 0: Strip CRLF line endings (do this first, always)

**Why:** FASTA files from many real sources (UniProt exports, anything
touched on Windows — including this repo's own original test database)
have Windows-style CRLF line endings. Left in, the stray `\r` characters
get indexed by crux as part of the actual protein sequences and end up
embedded inside tryptic peptides at line-wrap boundaries — this silently
corrupts every downstream tab-delimited output in a way that's invisible
until something tries to parse a field expecting a clean value. Confirmed
by hitting this for real while building this pipeline. CD-HIT does
**not** strip this for you — verified: it passes `\r` straight through.

```bash
tr -d '\r' < your_database.fasta > your_database_clean.fasta
```

**Check:**
```bash
grep -c $'\r' your_database_clean.fasta   # should print 0
```

---

## Step 1: CD-HIT — collapse near-identical sequences

**Why:** Duplicate/near-duplicate protein entries make the search engine
double-count/misreport shared peptides. Real result on `PSM_BL.fasta` (22
sequences): 2 collapsed at `-c 0.95`, 20 remain (at the stricter `-c 1.0`,
exact-matches-only, 0 get removed — this particular database has no
100%-identical entries, only two 95%+ near-identical ones).

```bash
cd-hit -c 0.95 -n 5 -d 0 -M 0 -T 0 \
    -i your_database_clean.fasta \
    -o your_database_dedup.fasta
```

**Parameters:**
- `-c 0.95` — identity threshold to cluster together: sequences that are
  ≥95% identical get merged into one. (CD-HIT's own documented default is
  `-c 0.9`, looser than this — worth knowing if you ever omit `-c` by
  accident. Raise it to `-c 1.0` for exact-duplicates-only; lower it for
  more aggressive collapsing, common with messy genome annotations.)
- `-n 5` — word length for CD-HIT's clustering *heuristic*, not the
  identity calculation itself. CD-HIT doesn't align every pair of
  sequences directly (too slow for large databases) — it first indexes
  every sequence by all its overlapping length-`n` substrings ("words":
  at `n=5`, every 5-amino-acid stretch), and only runs a real alignment on
  pairs that already share enough of those words to *possibly* meet the
  `-c` threshold. Two sequences that are ≥95% identical are mathematically
  guaranteed to share long exact-match stretches, so a 5-length word is
  short enough to always find that overlap; at a much lower `-c` (more
  differences allowed, so exact-match stretches get shorter and rarer),
  you need a *shorter* word (`-n` 2-4) or the pre-filter would wrongly
  throw out real matches before they're ever aligned. So `-n` doesn't
  change *which* sequences end up in a cluster (`-c` and the real
  identity computation decide that) — it only changes how CD-HIT
  efficiently narrows down candidate pairs to check. That's also why
  CD-HIT enforces the pairing (`cd-hit -h` for the full `-c`→`-n` table)
  and rejects a mismatched combination outright rather than silently
  giving you a wrong/incomplete clustering.
- `-d 0` — keep the full original header in the output (default truncates
  at the first space, which drops the protein description/`GN=` gene tag
  you'll likely want later).
- `-M 0 -T 0` — no memory cap, use all available CPU threads.

**Check:**
```bash
grep -c '^>' your_database_clean.fasta   # sequences going in
grep -c '^>' your_database_dedup.fasta   # sequences coming out

# Show only clusters that actually merged >1 sequence (most clusters are
# singletons and aren't interesting here)
awk '/^>Cluster/{if(block!="" && n>1) print block; block=$0"\n"; n=0; next}
     {block=block $0"\n"; n++}
     END{if(n>1) print block}' your_database_dedup.fasta.clstr
```
Real output on `PSM_BL.fasta`:
```
>Cluster 11
0	221aa, >sp|B5A7N4|XYN2_TRIHA... *
1	190aa, >sp|P48793|XYN_TRIHA... at 99.47%

>Cluster 12
0	147aa, >sp|P69892|HBG2_HUMAN... *
1	147aa, >sp|P69891|HBG1_HUMAN... at 99.32%
```
(the line ending in `*` is the cluster's representative sequence, kept in
the output; other lines were merged into it and dropped)

CD-HIT also prints a summary to stdout ending in `program completed !` —
if that line is missing, something failed; scroll up for the actual error.

---

## Step 2: Build the search database and index it

**Why:** Real samples almost always contain lab contaminants (keratins,
trypsin, BSA) regardless of what you're studying. Skipping this doesn't
mean those peptides go unidentified — they get **misassigned to whatever
real protein in your database looks similar**, a real source of wrong
identifications. `reference_data/common_contaminants.fasta` (116
sequences, the standard cRAP database) is bundled in this repo for
exactly this.

```bash
cat your_database_dedup.fasta "$CONTAM" > your_database_final.fasta

"$CRUX" tide-index \
    --decoy-format peptide-reverse \
    --overwrite T \
    --mods-spec C+57.02146 \
    --output-dir tide_index_log \
    your_database_final.fasta tide_index
```

**Parameters:**
- `--decoy-format peptide-reverse` — generates decoys by reversing each
  target peptide, for target-decoy FDR estimation in percolator later.
  (Older tutorials use `protein-reverse`; current crux 4.x renamed the
  options to `none`/`shuffle`/`peptide-reverse`.)
- `--mods-spec C+57.02146` — fixed carbamidomethylation on cysteine, the
  standard modification from iodoacetamide alkylation during sample prep.
  Drop/change this if your samples weren't alkylated.
- Decoys are generated automatically here — you do **not** need a
  separate decoy-database-building step (older workflows had you do this
  by hand with tools like OpenMS `DecoyDatabase`).

**Check:**
```bash
tail -20 tide_index_log/tide-index.log.txt
```
Look for:
- `Cleaved N protein sequences in total` — should equal your combined
  sequence count (`grep -c '^>' your_database_final.fasta`).
- `Ignoring N peptide sequences containing unrecognized characters` — a
  *small* N (e.g. 0-2) is normal, from real non-standard amino acid codes
  (`U` = selenocysteine, `B`/`Z` = ambiguity codes) that occasionally show
  up in real protein databases. A **large** N here means something's
  actually wrong with the input — check for leftover CRLF (step 0) first.
- `Return Code:0` at the end — anything else means it failed.

---

## Step 3: ThermoRawFileParser — RAW → mzML

```bash
mkdir -p mzML
"$TRFP" -i=your_run.raw -o=mzML -f=2 -l=2 2>&1 | tee mzML/your_run.convert.log
```

**Parameters:**
- `-f=2` — **indexed** mzML output. Always use this, not `-f=1` (plain
  mzML) — tools that look up specific spectra by ID (percolator's RT
  lookup, moFF if you go on to quantify) rely on the byte-offset index for
  fast access; without it, the same lookup can take 100x+ longer.
- `-l=2` — log verbosity (info level); needed to actually see per-spectrum
  warnings in the output.

**Check:**
```bash
tail -5 mzML/your_run.convert.log                       # did it finish, how many errors/warnings?
grep -c "has no m/z data" mzML/your_run.convert.log      # empty-spectrum count
grep -c '<spectrum id=' mzML/your_run.mzML               # total spectrum count
```
`RAW file cannot be processed since it is still being acquired` means that
specific `.raw` file was never cleanly closed by the instrument software
(a real flag inside the file, not a parser bug) — exclude that file and
move on, don't try to force past it.

A handful of "no m/z data" warnings among thousands of spectra is normal
(empty MS2 scans on low-signal precursors). If it's a large fraction,
check *which MS level* is affected before assuming something's broken —
a high empty-MS2 rate is usually still fine, but a high empty-**MS1**
rate is worth a closer look at *where in the run* (by retention time)
those empty scans fall: clustered at the very start/end of the run is
plausibly dead time (benign); a contiguous block in the middle, clean
before and after, points at a real, temporary acquisition event (e.g. an
electrospray dropout) rather than something to shrug off.

---

## Step 4: crux tide-search

```bash
mkdir -p search_out
"$CRUX" tide-search \
    --concat T \
    --precursor-window 50 \
    --mz-bin-width 1.0005079 \
    --overwrite T \
    --output-dir search_out \
    mzML/your_run.mzML tide_index
```

Pass multiple mzML files at once (space-separated, before the index path)
to search them together into one combined `tide-search.txt`.

**Parameters:**
- `--concat T` — write targets and decoys to one combined output file
  (required for percolator's `--search-input concatenated` in the next
  step).
- `--precursor-window 50` — ±50 ppm precursor mass tolerance. Widen for
  lower-resolution instruments, narrow for high-res Orbitrap/tribrid data
  if you know your instrument's real mass accuracy is tighter.
- `--mz-bin-width 1.0005079` — fragment-ion bin width used to discretize
  spectra before scoring. **Get this backwards and it's easy to do, so
  double-check:** crux's own `--help` text says *"For low resolution ion
  trap ms/ms data 1.0005079 and for high resolution ms/ms 0.02 is
  recommended. Default = 0.02."* — i.e. `1.0005079` is the **low-res**
  value, and crux's actual default is the **high-res** one (`0.02`), not
  `1.0005079`. `1.0005079` is used here because this particular test
  dataset really is low-resolution ion trap data (CID fragmentation, not
  HCD) — confirmed directly from the mzML, see the callout below, don't
  just assume. If your own data is high-resolution (Orbitrap/HCD), use
  `--mz-bin-width 0.02` instead — using the low-res bin width on
  high-resolution data throws away real discriminating power in the
  search. Crux can also estimate this for you automatically:
  `--auto-mz-bin-width warn` (falls back to the default if estimation
  fails) or `fail` (errors instead of silently guessing wrong).

**How to tell high- vs. low-resolution from an mzML file, without
guessing:** ThermoRawFileParser embeds the instrument's own per-scan
"filter string" (cvParam accession `MS:1000512`) — Thermo's own scan
description, straight from the instrument. `FTMS` = Orbitrap/FT
(high-res); `ITMS` = ion trap (low-res). Check across the whole file, not
just one scan, since some methods deliberately mix them (e.g. high-res
MS1 survey scans + fast low-res ion-trap MS2 for speed):
```bash
grep -o 'accession="MS:1000512" value="[^"]*"' your_run.mzML \
    | sed 's/.*value="//;s/"$//' | cut -c1-4 | sort | uniq -c
```
Real output on this repo's test data — every single scan, MS1 and MS2
alike, is `ITMS`:
```
    5574 ITMS
```
(a mixed-mode file would show both `FTMS` and `ITMS` counts; match your
`--mz-bin-width` to whichever one dominates your MS2 scans specifically,
since that's what's actually being bin-discretized for scoring)

**Check:**
```bash
grep -E "spectrum-charge combinations (loaded|searched)|Return Code" search_out/tide-search.log.txt
wc -l search_out/tide-search.txt
```
`N spectrum-charge combinations loaded, M spectrum-charge combinations
searched` — M is normally much smaller than N (most loaded
spectrum-charge combinations get filtered out before a full search, e.g.
as duplicates of another charge state or precursor-mass mismatches — this
is expected, not a sign of data loss).

---

## Step 5: percolator

```bash
"$CRUX" percolator \
    --test-fdr 0.1 \
    --train-fdr 0.1 \
    --search-input concatenated \
    --overwrite T \
    --output-dir search_out \
    search_out/tide-search.txt
```

**Parameters:**
- `--test-fdr` / `--train-fdr` — the FDR percolator uses to *train* its
  scoring model (10% here, matching this lab's exploratory/QC default).
  **Important:** this does NOT filter the output file (see check below) —
  it only affects model training.
- For a real study, 1% (`--test-fdr 0.01 --train-fdr 0.01`) is the
  standard reported in publications. You can also filter the existing
  output at a stricter threshold without rerunning percolator (see step
  6) — only rerun percolator if you want the model itself retrained at a
  different FDR.

**Check:**
```bash
tail -15 search_out/percolator.log.txt   # "Return Code:0" = success
wc -l search_out/percolator.target.psms.txt
```
**`percolator.target.psms.txt` contains every target PSM, not just
confident ones** — despite the FDR flags above. It has a `q-value` column
(4th column) you filter yourself in the next step. A file with hundreds or
thousands of rows can easily have zero rows that pass a real FDR
threshold — that's not a bug, check the actual confident count in step 6
before assuming something's wrong.

---

## Step 6: Protein count and gene list

**Why here, not step 5:** the FDR threshold that actually determines
"confident" is applied here, by filtering percolator's `q-value` column —
not automatically by percolator itself (see step 5's check).

```bash
# Confident PSMs at 10% FDR (column 4 = q-value; NR>1 skips the header)
awk -F'\t' 'NR>1 && $4<=0.1' search_out/percolator.target.psms.txt | wc -l

# Unique confident protein IDs -> protein count
awk -F'\t' 'NR>1 && $4<=0.1 {print $7}' search_out/percolator.target.psms.txt \
    | tr ',' '\n' | sort -u > confident_proteins.txt
wc -l confident_proteins.txt
```

**Optional: map to gene symbols.** UniProt-style FASTA headers carry a
`GN=GeneSymbol` tag (e.g. `>sp|P02769|ALBU_BOVIN ... GN=ALB ...`) — this
pulls that back out for each identified protein:

```bash
while read -r pid; do
    grep -m1 -F ">${pid}" your_database_final.fasta | grep -oP 'GN=\K\S+'
done < confident_proteins.txt | sort -u > gene_list.txt
wc -l gene_list.txt
```

**Real caveat, not hypothetical:** this only works for entries whose
header actually has a `GN=` tag. The bundled `common_contaminants.fasta`
does **not** carry `GN=` tags at all (it's a stripped-down contaminant
list, not a full UniProt export) — so any confidently-identified
contaminant will correctly show up in `confident_proteins.txt` but won't
produce a line in `gene_list.txt`. That's expected, not a bug in the
command — `wc -l gene_list.txt` will legitimately be smaller than
`wc -l confident_proteins.txt` whenever contaminants are among your hits.
If you need gene symbols for contaminant hits too, source a
fully-annotated contaminant FASTA instead (e.g. MaxQuant's, which does
include gene names).

**Check for a stricter FDR without rerunning percolator** (e.g. the 1%
publication standard):
```bash
awk -F'\t' 'NR>1 && $4<=0.01 {print $7}' search_out/percolator.target.psms.txt \
    | tr ',' '\n' | sort -u | wc -l
```

---

## Step 7: Spectral counting (protein-level abundance, no MS1 needed)

**Why:** Step 6 gets you a protein *list*. `crux spectral-counts` goes one
step further — a lightweight per-protein abundance estimate based purely
on how many PSMs/spectra map to each protein, no MS1 peak extraction
required. It's a real crux subcommand, not something hand-rolled with
awk, and it applies the FDR threshold itself rather than you filtering by
hand first. This is not the same as MS1-intensity-based quantification
(moFF, covered separately) — spectral counting is coarser, but it's fast
and needs nothing beyond what you already have from step 5.

```bash
"$CRUX" spectral-counts \
    --measure NSAF \
    --threshold 0.1 \
    --protein-database your_database_final.fasta \
    --overwrite T \
    --output-dir spectral_counts \
    search_out/percolator.target.psms.txt
```

**Parameters:**
- `--measure NSAF` — Normalized Spectral Abundance Factor: each protein's
  spectral count divided by its length, then normalized so all proteins'
  values sum to 1 — this makes long proteins (which naturally rack up more
  PSMs just from having more possible peptides) comparable to short ones.
  `--measure RAW` gives plain spectral counts instead (no length
  normalization, and doesn't need `--protein-database` at all) if you just
  want the simplest possible number. `dNSAF`, `SIN`, and `EMPAI` are other
  supported measures — see `crux spectral-counts` (no arguments) for the
  full list.
- `--threshold 0.1` — the q-value cutoff, applied internally (same idea as
  step 6's `--max-q-value`, but crux does the filtering for you here).
  Match this to whatever FDR you actually want to report at (e.g. `0.01`
  for the 1% publication standard).
- `--protein-database` — required for any length-normalized measure
  (`NSAF`/`dNSAF`/`SIN`/`EMPAI`); not needed for `RAW`. Use the same
  `..._final.fasta` you built in step 2 (dedup + contaminants) — it needs
  to match what you actually searched against.

**Check:**
```bash
tail -8 spectral_counts/spectral-counts.log.txt
cat spectral_counts/spectral-counts.target.txt
```
Real output on this repo's test data (single-file run):
```
INFO: Number of matches: 532
INFO: Number of matches passed the threshold: 10
INFO: Number of peptides: 10
INFO: Number of proteins: 10
INFO: Return Code:0
```
```
protein id	NSAF
sp|CYC_HORSE|	0.28302339
sp|NQO2_HUMAN|	0.12920633
sp|CAH2_BOVIN|	0.1142979
...
```
The `Number of proteins: N` line in the log is your final protein count —
cross-check it against `wc -l spectral_counts/spectral-counts.target.txt`
(subtract 1 for the header) as a sanity check; they should match. If they
don't, or `Return Code` isn't `0`, something's wrong — check the log for
the actual error rather than trusting the output file blindly.

(Optional: `--parsimony simple` or `--parsimony greedy` collapses proteins
whose identified peptides are a subset of/shared with a longer protein
into one group, and adds a `parsimony rank` column — worth knowing about
for real studies with a large database where the same peptide often maps
to multiple homologous proteins; skip it for a small demo database like
this one, where it won't do much.)

---

## Full quick-reference (all commands, no explanation)

```bash
conda activate proteomics

tr -d '\r' < your_database.fasta > your_database_clean.fasta

cd-hit -c 0.95 -n 5 -d 0 -M 0 -T 0 -i your_database_clean.fasta -o your_database_dedup.fasta

cat your_database_dedup.fasta "$CONTAM" > your_database_final.fasta
"$CRUX" tide-index --decoy-format peptide-reverse --overwrite T --mods-spec C+57.02146 \
    --output-dir tide_index_log your_database_final.fasta tide_index

mkdir -p mzML
"$TRFP" -i=your_run.raw -o=mzML -f=2 -l=2 2>&1 | tee mzML/your_run.convert.log

mkdir -p search_out
"$CRUX" tide-search --concat T --precursor-window 50 --mz-bin-width 1.0005079 \
    --overwrite T --output-dir search_out mzML/your_run.mzML tide_index

"$CRUX" percolator --test-fdr 0.1 --train-fdr 0.1 --search-input concatenated \
    --overwrite T --output-dir search_out search_out/tide-search.txt

awk -F'\t' 'NR>1 && $4<=0.1 {print $7}' search_out/percolator.target.psms.txt \
    | tr ',' '\n' | sort -u > confident_proteins.txt
wc -l confident_proteins.txt

"$CRUX" spectral-counts --measure NSAF --threshold 0.1 \
    --protein-database your_database_final.fasta \
    --overwrite T --output-dir spectral_counts \
    search_out/percolator.target.psms.txt
```
