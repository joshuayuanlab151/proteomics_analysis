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
   database is intentionally tiny — 22 proteins — so don't expect large
   numbers of hits from it).
3. A modification setting mismatch, e.g. `--mods-spec C+57.02146` assumes
   cysteine alkylation was performed; if it wasn't, remove that flag when
   rebuilding the index.
4. The RAW→mzML conversion silently produced near-empty files — check file
   sizes (doc 3.3) before assuming the search itself is broken.
