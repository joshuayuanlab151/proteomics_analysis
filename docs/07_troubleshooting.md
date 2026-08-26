# 7. Troubleshooting

Real problems hit while building and testing this pipeline, and their
fixes — read this before assuming something is broken beyond repair.

*(Need more compute than your own laptop has for a real study? See
[doc 8](08_hellbender_slurm.md) for running this on Mizzou Hellbender.)*

## "`pip install moFF` doesn't work" / "moFF does something completely unrelated to proteomics"

**Don't use `pip install moFF`.** PyPI hosts a package named `moFF` (and
`MOFF`) that predicts CRISPR/Cas9 off-target effects — a different tool
from a different group, unrelated to mass spectrometry, that happens to
share a name. `import moff` from that package will succeed and give you no
error, just wrong results/nonsense if you tried to use it as if it were the
proteomics tool. Always install the real moFF via
`00_setup_environment.sh` (in `mac_linux/scripts/` or `windows/scripts/`),
which clones
[compomics/moFF](https://github.com/compomics/moFF) from GitHub and runs
its scripts (`moff_all.py`) directly — it is not designed to be pip-
installed as a package at all.

## `crux percolator` crashes with "Segmentation fault"

This pipeline downloads crux directly rather than through conda (see doc
1.2) and pins it to **crux-4.3** (2025-03-25), not the newer **4.3.2**
release. That's a specific, tested finding, not general version-caution:
4.3.2's bundled percolator reliably segfaulted partway through training on
this pipeline's example data, while 4.3 — the release immediately before
it — completed the identical run cleanly and produced identical, correct
output all the way through moFF. So this is an isolated regression in
4.3.2 specifically, not a reason to distrust the 4.3.x line in general.
crux-4.2 was also tested and works, but 4.3 is preferred since it's newer
and (as a bonus) includes retention time directly in its output.

If you ever change `CRUX_VERSION` in `00_setup_environment.sh` (in either
`mac_linux/scripts/` or `windows/scripts/` — update both, they're
independent copies) to a different release (e.g. once a fixed 4.3.3+ comes
out), re-run the full pipeline on the example data first and confirm
`percolator.target.psms.txt` actually gets produced before trusting it on
real data — don't assume a newer version number is automatically safe.

## `crux: /lib64/libm.so.6: version 'GLIBC_2.29' not found` (common on Hellbender / older HPC Linux)

The `mac_linux/` track downloads crux as a prebuilt binary from GitHub
rather than through conda. That binary is compiled against a fairly recent
glibc/libstdc++. Some Linux systems — including HPC clusters running
older, deliberately-conservative OS images, such as Mizzou Hellbender —
ship an older glibc than that, so the binary refuses to run and prints
something like:

```
crux: /lib64/libm.so.6: version `GLIBC_2.29' not found (required by crux)
crux: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.26' not found (required by crux)
```

`00_setup_environment.sh` checks for exactly this failure after
downloading crux and works through fallback tiers on Linux, in order.

*(Tier numbers below match what the script actually prints, e.g. "Tier 2:
falling back to..." — tier 1 is just the plain downloaded binary that got
us into this mess.)*

**Tier 2 (automatic): install crux from bioconda, pinned to the 4.1
build** (`conda create -n crux-fallback -c bioconda
crux-toolkit=4.1=h503566f_3`) — this gets crux running immediately and
automatically, no manual steps needed. **It is not the recommended
end state, though:** 4.1 is two releases behind the current 4.3 and has
known bugs of its own (see the CV-term entry below, which this tier's
wrapper works around by sanitizing mzML files rather than fixing at the
source). **If you can spare the 30-90 minutes, prefer Tier 4 (source
build) below, which gets you the actual current release (4.3)** — best
practice is to run the newest crux version you can, not settle for
whatever came up first. Tier 2 exists so you're not blocked while you
decide whether to invest that time, and remains a reasonable choice for
a quick one-off run.

The pinned build number matters here: bioconda's *newest*
crux-toolkit build (4.2) has the exact same problem as the direct
download — it happened to get compiled on a newer CI base image, even
though 4.1 and 4.2 are sibling releases (tagged 8 minutes apart, same
day, identical CLI). We checked by downloading each bioconda build and
inspecting the actual versioned symbols the `crux` binary inside it
requires (`strings crux | grep -oE 'GLIBC(XX)?_[0-9.]+' | sort -V`):

| build | max GLIBC needed | max GLIBCXX needed |
|---|---|---|
| 4.2 (`h9ee0642_0`) | 2.29 | 3.4.26 — same failure as the direct download |
| 4.1 (`h503566f_3`, newest 4.1 build) | 2.14 | 3.4.19 — comfortably under AlmaLinux 8.9's glibc 2.28 |

So pinning to `crux-toolkit=4.1=h503566f_3` specifically (not just "4.1",
in case a future 4.1 rebuild changes this) is what makes tier 2 actually
work on Hellbender, rather than failing the same way tier 1 did.

If you'd already hit this before this fix landed, `00_setup_environment.sh`
may have left behind a `crux-fallback` conda environment with the broken
4.2 build in it. The script now checks the installed version and rebuilds
that environment automatically if it's not 4.1 — you don't need to
manually `conda env remove -n crux-fallback` first, just re-run
`./00_setup_environment.sh`.

**Tier 3 (automatic, rarely needed now): run crux inside a container**
via Apptainer/Singularity, pulling
`docker://quay.io/biocontainers/crux-toolkit:4.2--h9ee0642_0`. Only
reached if tier 2 somehow still fails. A container bundles its own
glibc/libstdc++ inside its own filesystem, so it doesn't matter what the
host OS provides — but it does need Singularity/Apptainer to be available
to your account, which is not guaranteed (Hellbender documents "limited
support for Singularity for advanced users"). The script looks for
`singularity` or `apptainer` on `PATH`, and if neither is there, tries
`module load singularity` / `module load apptainer` before giving up. If
this tier fails because neither command/module can be found, email
`itrss-support@umsystem.edu` and ask them to confirm Singularity/Apptainer
is enabled on your account.

**Tier 4 (manual, recommended over Tier 2 when you have the time): build
crux from source.** This is how you get the current release (4.3) working
on Hellbender — the only tier that does. Compiling crux
against the system's own compiler and libraries sidesteps the glibc
mismatch entirely, by construction. This is a real C++ build, though, not
"download and unzip": it needs `cmake`, a C++ compiler, and system
development headers (openssl, zlib, libcurl, libuuid, libpulse). On a
shared HPC system, some of these may only exist as versioned environment
modules (e.g. `cmake/3.26.3_gcc_9.5.0`, not bare `cmake`) rather than on
`PATH` by default.

**Despite what crux's own generic `INSTALL` file says, `svn`/subversion
is NOT actually required.** Confirmed by tracing the entire default build
chain (top-level `CMakeLists.txt` → `ext/CMakeLists.txt` →
`download_pwiz.cmake` → `build_pwiz.cmake`) end to end: zero
svn/subversion references anywhere in it. ProteoWizard is fetched via a
plain HTTPS URL (a pre-packaged source tarball), not an SVN checkout —
the `INSTALL` file's mention of svn is almost certainly a stale holdover
from when crux's own source lived on SourceForge SVN, not an actual
build-time dependency of the current CMake-based build.
`build_crux_from_source.sh` no longer checks for it. If you were
previously told (by this doc, in an earlier version) that missing `svn`
blocks this tier — it doesn't; try again.

It can also take 30-90+ minutes with continuous network access, so run it
as its own multi-hour SLURM job — **not** inside the short interactive
session used for normal setup:
```bash
cd $HOME/data/proteomics_pipeline_course/mac_linux/scripts/slurm
mkdir -p logs
sbatch build_crux_from_source.slurm
```
(`mac_linux/scripts/build_crux_from_source.sh` is the underlying script,
in case you want to run it directly in an interactive session with enough
time instead.)

Whichever tier ends up working, `software/crux/bin/crux` ends up being
either a real compiled binary (tier 4) or a small wrapper script that
activates the conda env or execs into the container (tiers 2/3) — either
way, every other script (`02_build_database_index.sh`,
`03_tide_search.sh`, `04_percolator.sh`, ...) keeps working completely
unmodified, since they all just call `software/crux/bin/crux` and don't
know or care which crux is actually behind it. The tier 2 wrapper also
loads the `miniconda3` module itself if `conda` isn't already on `PATH`
in whatever shell calls it, so you never need to manually
`module load miniconda3` / `conda activate crux-fallback` before running
`software/crux/bin/crux` (or any pipeline script) yourself — just call it
directly, same as if it were the real binary. `00_setup_environment.sh`'s
final "Setup complete" message also reports which backend is actually
running, by inspecting `software/crux/bin/crux` itself.

Tiers 2 and 3 end up using crux **4.1/4.2** instead of the current
**4.3** — only Tier 4's source build gets you 4.3 on Hellbender. For this
pipeline's own moFF-bound output, the version gap happens to be low-risk:
`05_percolator_to_moff.py` computes retention time itself directly from
the mzML files rather than relying on tide-search's 4.3-only
`retention time` column, so that specific output is identical either way.
That's a reason not to panic if you're stuck on Tier 2 for now, not a
reason to prefer it — general bug fixes between 4.1 and 4.3 (like the
CV-term issue above) still only get fixed by actually running 4.3, so
Tier 4 remains the recommended target when you have the time for it.

If setup was interrupted partway through, or you want to force it to
retry, just delete `software/crux/bin/crux` and re-run
`00_setup_environment.sh` — it re-downloads and walks the tiers again
(reusing the `crux-fallback` conda env if it's already correct, and any
already-pulled `crux.sif`, so the retry is quick), rather than getting
stuck on whatever failed last time.

If you ever see this and want to confirm which tier actually kicked in,
check `software/crux/bin/crux` — a real downloaded/built binary is
several MB, a fallback wrapper is a tiny few-line shell script, and
`conda list -n crux-fallback crux-toolkit` (if tier 2) will show the
exact installed version.

## `FATAL: An exception occurred: Invalid cvParam accession "1003096"` (or similar) from `crux tide-search`

A direct consequence of the tier 2 (conda) fallback above, not a separate
bug: the pinned bioconda build is **crux 4.1, from October 2021**. Its
bundled mzML parser only recognizes CV (controlled vocabulary) terms that
existed in the PSI-MS ontology as of that date. Current
ThermoRawFileParser writes a couple of terms that didn't exist yet in
2021, and crux 4.1 hard-fails on any accession it doesn't recognize
rather than skipping it.

**Confirmed for real, exhaustively, on actual Hellbender data:** exactly
two such terms show up, and only two — `MS:1003096` (the specific
instrument model, e.g. "LTQ Orbitrap Velos Pro") and `MS:1003145`
(ThermoRawFileParser identifying itself as the conversion software).
Both are one-off, file-level metadata that appears exactly once per file
(inside `referenceableParamGroupList` and `softwareList` respectively) —
**not** something that recurs per-spectrum, so this isn't an open-ended
problem despite how alarming a bare "invalid accession" error looks.
Checked directly (`grep -o 'accession="MS:100[3-9][0-9]*"' yourfile.mzML
| sort -u`) against every peak/spectrum-level field (the ones appearing
thousands of times — m/z arrays, filter strings, scan times, etc.) and
they're all long-standing `1000xxx`-range terms, comfortably older than
crux 4.1.

**`03_tide_search.sh` handles this automatically now:** it detects the
tier 2 conda-fallback crux (same detection used elsewhere: checking for
`crux-fallback` in the wrapper script) and, only in that case, searches
*sanitized copies* of your mzML files instead of the originals — your
real converted mzML files are never modified. The two known terms are
patched: `MS:1003096` is removed outright (pure descriptive metadata,
safe to drop), and `MS:1003145` is replaced with the long-standing,
generic `MS:1000799` ("custom unreleased software tool", whose own
definition literally says "the value should describe the software") with
`value="ThermoRawFileParser"` — swapped rather than deleted, since it's
the only cvParam in its `<software>` element and is referenced elsewhere
via `softwareRef`, so deleting it risked trading this error for a
different structural one. This has **not** been verified against a real
run of crux 4.1 (no crux 4.1 available to test outside Hellbender itself)
— it's verified to be well-formed XML with the correct structure and a
semantically appropriate substitute term, and crux tide-search doesn't
use the mzML's own internal `<indexList>` (it converts to its own
internal "spectrumrecords" format), so editing header bytes ahead of that
index is harmless for this specific use.

If you still hit an "Invalid cvParam accession" error after this fix (a
different accession number), that means a third term exists that wasn't
in this dataset's files — report the exact number and where in the file
it appears (`grep -n 'accession="MS:<the number>"' yourfile.mzML`) so it
can be patched precisely rather than guessed at, the same way the first
two were found and confirmed.

If you're not on tier 2 (i.e. `crux -v` shows something other than
`4.1-...-2021...`), this doesn't apply to you at all — this is
specifically a byproduct of the old crux 4.1 build, not a general crux or
ThermoRawFileParser problem.

## A `crux-output` folder keeps appearing in `scripts/`

Harmless. crux always creates a default `crux-output` log directory in
whatever folder you ran it from, regardless of `--output-dir` — this
appears to be unconditional crux behavior, not something our scripts
control. `02_build_database_index.sh`, `03_tide_search.sh`, and
`04_percolator.sh` all clean it up automatically after each run; if you
ever invoke crux directly by hand instead, you may see it appear and can
safely delete it.

## Installing the `moff` conda environment seems to hang forever

moFF pins old, exact package versions (Python 3.6, numpy 1.19.5, pandas
1.1.5, etc.). Asking conda's default solver to satisfy all of those
constraints simultaneously against modern channel indexes can take an
extremely long time (sometimes 15+ minutes, sometimes effectively forever)
— this is a known conda pain point with old pinned dependency sets, not a
problem with your setup.

`00_setup_environment.sh` avoids this by creating a plain `python=3.6`
environment (a fast, trivial solve) and then `pip install`-ing the exact
pinned versions directly — pip's resolver handles this case much faster
than conda's. If you're troubleshooting this by hand, do the same:

```bash
conda create -n moff python=3.6
conda activate moff
pip install numpy==1.19.5 pandas==1.1.5 scipy==1.5.3 scikit-learn==0.24.1 \
    pymzml==2.4.7 pyteomics==4.4.2 simplejson==3.17.2 \
    brain-isotopic-distribution==1.5.3 lxml
```

(`pynumpress` is in moFF's official environment file but is optional — see
next item — and its build can fail on some machines, so it's intentionally
left out above.)

## `pynumpress` fails to build ("xcrun: error... missing compatible architecture")

This happens on Apple Silicon Macs because `pynumpress` needs to compile a
C extension and the toolchain available can't cross-compile for the
x86_64 conda environment. **`pynumpress` is optional** — it only matters
for reading numpress-compressed mzML files, which ThermoRawFileParser does
not produce by default. Just skip it; moFF runs fine without it.

## `No module named 'brainpy._c.composition'` warning when running moFF

Harmless. `brain-isotopic-distribution` has an optional compiled
accelerator; without it, it falls back to a pure-Python implementation.
moFF still works correctly, just marginally slower.

## `crux tide-index` errors: `Invalid value for 'decoy-format'`

Older tutorials/documentation reference `--decoy-format protein-reverse`.
Current crux (4.x) renamed the available options to `none`, `shuffle`, and
`peptide-reverse`. Use `--decoy-format peptide-reverse` (already the
default in `02_build_database_index.sh`).

## ThermoRawFileParser conversion gets killed instantly (exit code 137) on macOS

This applies to the `mac_linux/` track on macOS specifically (not Linux,
not Windows). macOS refuses to execute a freshly downloaded, unsigned
ARM64/x64 binary outright (silently killing it, even without a quarantine
flag set) unless it's been code-signed. Fix by ad-hoc signing it once,
from inside `mac_linux/`:

```bash
codesign --force --deep -s - software/ThermoRawFileParser/ThermoRawFileParser \
    software/ThermoRawFileParser/*.dylib
codesign --force --deep -s - software/crux/bin/crux
```

`mac_linux/scripts/00_setup_environment.sh` already does this
automatically — you should only need this if you moved/redownloaded a
binary manually. Not needed on Linux or Windows.

## `RAW file cannot be processed since it is still being acquired`

This is a real flag stored inside that specific `.raw` file (`InAcquisition`
in the Thermo file format), not a parser bug — confirmed by testing: the
same tool, same machine, and same command convert every other file from
the same batch without issue. It means the instrument software never
marked that particular acquisition as cleanly finished, typically because
the run was stopped early/aborted, or the file was copied off the
instrument before Xcalibur finished writing/closing it.

**`01_convert_raw_to_mzml.sh` skips a file that fails like this and keeps
converting the rest of the batch** — it doesn't stop the whole run over
one flagged file. At the end it prints a summary of any files that failed
and were skipped, and each file's full output is saved to
`<mzml_dir>/<name>.convert.log` for reference. The rest of the pipeline
then runs normally on whichever files did convert successfully.

There is no supported command-line flag to force past this check (it's not
gated by `--ignoreInstrumentErrors`). Options:

- Ask whoever ran the acquisition whether that run was intentionally
  stopped early — if so, the file may simply need to be excluded from
  analysis.
- Check if a properly closed copy of the same file exists elsewhere (e.g.
  it was already converted once via vendor software or a different
  ThermoRawFileParser install at acquisition time).
- If you have access to Proteome Discoverer / Xcalibur on the acquiring
  instrument's Windows PC, re-exporting a copy from there sometimes
  resolves the flag.

## `WARN Spectrum #N has no m/z data` during RAW→mzML conversion

Traced directly to ThermoRawFileParser's source
(`Writer/MzMlSpectrumWriter.cs`): it logs this whenever a specific scan's
peak array comes back empty (`mzData.masses.Length == 0`) — i.e. spectrum
`#N` in your RAW file genuinely has zero detected peaks. This is normal
instrument behavior, not file corruption or a parser bug: it usually
means a data-dependent MS2 scan got triggered on a precursor that turned
out to have no significant fragment ion signal above threshold. The
spectrum still gets written to the mzML (as an empty spectrum) rather
than the conversion failing.

**A handful of these among thousands of spectra is normal** — those scans
just won't identify anything in tide-search, same as any other low-quality
spectrum, and cost you nothing else. Only worth investigating if a large
fraction of spectra warn like this, which points at a different, real
problem instead.

`01_convert_raw_to_mzml.sh` now does this check for you automatically:
it saves each conversion's output to `<mzml_dir>/<name>.convert.log`
alongside the mzML file, and prints a one-line summary (count of empty
spectra out of the total) right after each file finishes converting if
any showed up. To check by hand, e.g. for an older log or a manual
conversion:
```bash
grep -c "has no m/z data" your_output_dir/yourfile.convert.log
grep -c '<spectrum id=' your_output_dir/yourfile.mzML
```
Checking the mzML directly (`defaultArrayLength="0"` is the exact
condition that triggers the warning) works even without a saved log:
```bash
grep -c 'defaultArrayLength="0"' yourfile.mzML
```
If most spectra are affected, that's more likely the "still being
acquired" / near-empty file issue above, not this one — check file size
and total spectrum count before assuming the RAW file itself is fine.

### It's more than "a handful" (e.g. 5-10%+) — is that actually a problem?

The total count alone can't tell you that — what matters is **which MS
level** the empty spectra are.

`check_mzml_ms_levels.py` (in `mac_linux/scripts/` and `windows/scripts/`,
identical) breaks the count down by MS level so you can tell which case
you're in:
```bash
conda activate proteomics
python3 check_mzml_ms_levels.py yourfile.mzML
```
```
yourfile.mzML:
  MS1: 0/1203 empty (0.0%)
  MS2: 1011/9075 empty (11.1%)
```

**Empty MS2 at this kind of rate is normal.** Aggressive data-dependent
acquisition (e.g. Top20/Top30 methods) routinely triggers fragmentation
on precursors that turn out to be near the noise floor, and
low-abundance/complex real samples can easily produce a 5-15% empty-MS2
rate without anything being wrong.

**Empty MS1 needs a closer look than the flat percentage gives you.** We'd
originally written that MS1 scans should essentially never be empty and a
high rate always means a real acquisition problem, then found real
counter-examples in this pipeline's own original 12-file test set (several
files with 18-77% empty MS1, clean MS2) that seemed to contradict that.
But looking at *where in the run* those empty scans actually fall settled
it — `check_mzml_ms_levels.py` now reports this automatically whenever
MS1 empty rate is >=5%:
```
  MS1: 3459/5039 empty (68.5%)  <-- worth a closer look, see retention-time breakdown below
  MS1 empty rate by retention time (2-min bins):
       0.0-  2.0 min:    0/  29 empty (  0.0%)
       ...
       6.0-  8.0 min:  239/ 261 empty ( 91.6%) ##################
       8.0- 10.0 min:  658/ 669 empty ( 98.4%) ###################
      10.0- 12.0 min:  656/ 670 empty ( 97.9%) ###################
      ...
      16.0- 18.0 min:  587/ 635 empty ( 92.4%) ##################
      18.0- 20.0 min:    0/  30 empty (  0.0%)
       ... (clean for the rest of a 110-minute run)
```
Two genuinely different patterns, and they mean different things:

- **Empty MS1 clustered at the very start and/or end of the run** —
  plausibly dead time (column equilibration before the gradient moves
  analyte, wash/re-equilibration after) — benign.
- **A sharp, contiguous block of empty MS1 in the MIDDLE of the run,
  clean before it and clean after it** (the real example above) — this is
  **not** dead time. Two different files in this pipeline's own test data
  showed almost exactly this pattern, in a similar retention-time window
  (~6-24 min), both bounded by long stretches of completely normal signal
  on either side. That points at a real, temporary acquisition event
  during that specific window — most plausibly an electrospray/ion-source
  dropout (clogged needle, air bubble, voltage instability) that
  self-resolved — not something to shrug off. Practically: anything that
  should have eluted during that window is essentially invisible for MS1
  quantification in that specific run. Worth flagging to whoever ran the
  acquisition, and worth deciding whether that run needs to be excluded or
  re-acquired if completeness in that RT range matters for your study.

Either way, this is diagnostic information for you to act on, not
something the pipeline itself filters or excludes automatically — a file
with a bad window like this still gets searched/quantified like any
other; nothing downstream treats it specially.

## moFF reports "Peak not detected" / 0 retrieved peaks for every single peptide

This was a **real bug in moFF itself**, not a config or data problem, found
and fixed while building this pipeline — worth understanding if you ever
see it come back (e.g. after a manual reinstall that skips the patch).

moFF's internal `scan_mzml()` / `pyMZML_xic_out()` functions used to fetch
each spectrum via pymzml's random-access `run[spectrum_id]` lookup. With
current pymzml + current ThermoRawFileParser-generated mzML, that call
reliably raises `AttributeError: 'NoneType' object has no attribute
'calling_instance'` for **every** spectrum — silently swallowed by a bare
`except: pass` in moFF's own code (which even has a code comment admitting
this random-access behavior is "a mystery"). The practical effect: moFF's
internal retention-time index ends up completely empty, so its search
window for every peptide is empty too, and it reports every single peptide
as "Peak not detected" — while running suspiciously fast, because it isn't
actually searching anything.

**`00_setup_environment.sh` fixes this automatically** (both tracks):
after cloning moFF, it applies `patches/moff_random_access_fix.patch`,
which replaces the broken random-access lookups with a single sequential
pass over the mzML file (building an in-memory `{spectrum_id: peaks}`
cache once per run instead of re-fetching per peptide) — same result, and
actually **faster**, not slower, than the original code.

If you ever reinstall moFF by hand instead of via the setup script, make
sure to also run (from inside `mac_linux/` or `windows/`):
```bash
cd software/moFF
git apply ../../../patches/moff_random_access_fix.patch
```
(`software/moFF` is two levels under `mac_linux/`, which is itself one
level under the repo root where `patches/` lives — three `../` in total.)
You'll know the patch didn't get applied if every peptide in
`*_moff.log` says "Peak not detected" / "Xic shape 0" despite your
identifications clearly being real (check with doc 5.4's sanity checks).

## Converting to indexed vs. plain mzML matters a lot for speed

Always convert with `-f=2` (indexed mzML), as `01_convert_raw_to_mzml.sh`
does by default — **not** `-f=1` (plain mzML). Tools that look up specific
spectra by ID (including moFF, even after the fix above) rely on the
byte-offset `<indexList>` at the end of an indexed mzML file for fast
lookups. Without it, the same lookup pattern falls back to scanning the
file from the start every time, which measurably turned a several-second
moFF run into a 25+ minute one on identical data during testing. If you
have existing plain (non-indexed) mzML files and moFF seems unreasonably
slow, that's the first thing to check — reconvert with `-f=2`.

## moFF quantification is slow even with everything above correct

Some slowness is still expected — moFF is pure Python and reads real MS1
signal for every identified peptide. If you piped its output through
something like `| tail`, remember `tail` (without `-f`) buffers until the
input stream ends, so you may see **zero output until the whole run
finishes** even though it's working the entire time — a shell buffering
quirk, not a hang. Run it directly in your terminal (not piped) so you can
see progress as it happens, and let it run in the background for larger
studies rather than waiting on it actively. Native Linux is also
meaningfully faster than an Apple Silicon Mac running this under Rosetta
emulation, if you're comparing speed across machines.

**Match-between-runs (`--mbr on`) is somewhat heavier than MBR off**, but
on the bundled 2-run example dataset the full pipeline (RAW conversion
through final quantification table) typically completes in well under 10
minutes on a normal, uncontended machine. It *can* run much slower (one
test run took over 15 minutes) if the machine is under heavy load from
other processes at the same time — moFF's MBR step is CPU-bound and
doesn't compete gracefully for cores. If a run seems unexpectedly slow,
check what else is using CPU before assuming something is broken. Real,
larger studies will take longer than the example regardless.

## My percolator results have very few / zero confident PSMs

Usually one of:
1. Wrong or mismatched database (wrong organism, or FASTA doesn't
   correspond to your sample).
2. Too strict an FDR for a very small test/database (the bundled example
   database is intentionally tiny — 22 real proteins, 138 once the
   required contaminant sequences are added, see doc 2.4 — so don't
   expect large numbers of hits from it).
3. A modification setting mismatch, e.g. `--mods-spec C+57.02146` assumes
   cysteine alkylation was performed; if it wasn't, remove that flag when
   rebuilding the index.
4. The RAW→mzML conversion silently produced near-empty files — check file
   sizes (doc 3.3) before assuming the search itself is broken.

## `percolator.target.psms.txt` has thousands of rows — is that all "confident"?

**No.** This tripped us up too: despite `04_percolator.sh` setting
`--test-fdr`/`--train-fdr`, `percolator.target.psms.txt` is **not**
filtered down to confident matches — percolator uses those FDR values to
*train* its scoring model, but still writes out every single target PSM
with a `q-value` column, whatever that value ends up being. A file with
thousands of rows can easily have zero rows that actually pass any
reasonable FDR threshold.

The actual filtering happens in the next step, `05_percolator_to_moff.py`'s
`--max-q-value` (default `0.1`, matching `04_percolator.sh`'s default) —
see docs 5.1/5.3/6.1. Before this was added, the pipeline was silently
passing every unfiltered PSM through to moFF for quantification, which
defeats the entire point of running percolator. Check
`Total PSMs written: X, skipped: Y (+ Z filtered out by q-value > ...)`
in that script's output to see how many actually survived filtering.

## `TypeError: float() argument must be a string or a number, not 'NoneType'` from `05_percolator_to_moff.py`

Traced to a real, previously-latent bug found while testing against real
data: **if your input FASTA (or the database prep step's contaminants
file) has Windows-style CRLF line endings**, the stray `\r` characters
get indexed by crux as part of the actual protein sequences and end up
embedded inside tryptic peptides at line-wrap boundaries. Downstream,
that `\r` lands inside a field of `percolator.target.psms.txt` (typically
inside a peptide sequence, right after the trailing `.` before the next
tab). Python's `csv` module treats a lone `\r` as a row terminator, so it
silently splits that one real row into two bogus ones — a garbled
short row with mostly-empty fields, which then breaks `float(row["q-value"])`.

You'd also have seen a symptom earlier, at the indexing step:
`INFO: Ignoring N peptide sequences containing unrecognized characters`
from `crux tide-index` — that N is exactly the CRLF damage, and is your
earliest warning something's wrong with the input FASTA's line endings.

**`02_build_database_index.sh` fixes this automatically** (both tracks,
as of this fix): it strips `\r` from both your input FASTA and the
bundled contaminants FASTA (`tr -d '\r'`) before doing anything else with
them, so this shouldn't happen anymore even if your source FASTA has CRLF
endings — which is common (UniProt exports, Windows-authored files; this
pipeline's own original test database had this exact issue). If you ever
reinstall/rebuild a database index by hand instead of via the script,
remember to strip `\r` yourself first:
```bash
tr -d '\r' < your_database.fasta > your_database_clean.fasta
```
