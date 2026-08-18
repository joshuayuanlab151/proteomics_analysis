# Linux command-line walkthrough: RAW → final protein/gene list

A back-to-back command reference for demoing the identification pipeline
by hand on Linux: CD-HIT → build the search database → ThermoRawFileParser
→ crux tide-search → percolator → protein/gene list. Every command and log
excerpt below was actually run against this repo's real test data
(`PSM_BL.fasta`, one real `.raw` file) while writing this doc — not
guessed.

This is a manual, step-by-step version of what `mac_linux/scripts/`
automates for you. Use this for explaining/demoing the mechanics; use the
scripts in `mac_linux/scripts/` for actually running real studies (they
also handle exact-duplicate removal without needing CD-HIT, for
Windows-track parity — see `docs/02_prepare_database.md`).

## Prerequisites

Assumes crux and ThermoRawFileParser are already installed (via
`mac_linux/scripts/00_setup_environment.sh`, or manually — see
`docs/01_install_software.md`), and the `proteomics` conda environment
exists. CD-HIT is not part of that environment by default — install it
once:

```bash
conda activate proteomics
conda install -y -c bioconda -c conda-forge cd-hit
```

Set these once per session so every command below can just be copy-pasted:

```bash
export CRUX=/path/to/mac_linux/software/crux/bin/crux
export TRFP=/path/to/mac_linux/software/ThermoRawFileParser/ThermoRawFileParser
export CONTAM=/path/to/reference_data/common_contaminants.fasta
export WORKDIR=/path/to/demo_run
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
by hitting this for real while building this pipeline (see
`docs/07_troubleshooting.md`). CD-HIT does **not** strip this for you —
verified: it passes `\r` straight through.

```bash
tr -d '\r' < your_database.fasta > your_database_clean.fasta
```

**Check:**
```bash
grep -c $'\r' your_database_clean.fasta   # should print 0
```

---

## Step 1: CD-HIT — collapse near-duplicate sequences

**Why:** Near-duplicate protein entries (common in some genome
annotations) make the search engine double-count/misreport shared
peptides. Real result on `PSM_BL.fasta` (22 sequences): 2 collapsed as
95%-identical, 20 remain.

```bash
cd-hit -c 0.95 -n 5 -d 0 -M 0 -T 0 \
    -i your_database_clean.fasta \
    -o your_database_dedup.fasta
```

**Parameters:**
- `-c 0.95` — identity threshold to cluster together (95%). Use `-c 1.0`
  for exact-duplicates-only; lower (e.g. `0.9`) to collapse more
  aggressively for very redundant genome annotations.
- `-n 5` — word size CD-HIT requires for `-c` in the 0.9-0.95 range (CD-HIT
  ties word size to identity threshold; see `cd-hit -h` for the mapping —
  wrong `-n` for a given `-c` is rejected outright, not silently wrong).
- `-d 0` — keep the full original header in the output (default truncates
  at the first space, which drops the protein description/`GN=` gene tag
  you'll likely want later).
- `-M 0 -T 0` — no memory cap, use all available CPU threads.

**Check:**
```bash
grep -c '^>' your_database_clean.fasta   # sequences going in
grep -c '^>' your_database_dedup.fasta   # sequences coming out
head -20 your_database_dedup.fasta.clstr # which sequences got merged, and with what
```
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
tail -20 tide_index_log/tide-index.log
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
see `mac_linux/scripts/check_mzml_ms_levels.py` and
`docs/07_troubleshooting.md` for why a flat percentage alone can be
misleading (real example there: a sharp mid-run dropout, not "always
bad").

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
- `--mz-bin-width 1.0005079` — standard fragment-ion bin width tuned to
  the mass of a proton; the conventional default for high-res MS2, rarely
  needs changing.

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
```
