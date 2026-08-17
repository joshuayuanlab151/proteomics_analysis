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
crux-toolkit=4.1=h503566f_3`) — this is what actually fixes it on
Hellbender. The pinned build number matters here: bioconda's *newest*
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

**Tier 4 (manual, last resort): build crux from source.** Compiling crux
against the system's own compiler and libraries sidesteps the glibc
mismatch entirely, by construction. This is a real C++ build, though, not
"download and unzip": it needs `cmake`, `svn` (crux's default build
fetches ProteoWizard as a tarball, not SVN — see the CMake logic in
`ext/CMakeLists.txt` — but `INSTALL` still lists both `cmake` and `svn` as
prerequisites, so treat both as required until proven otherwise), a C++
compiler, and system development headers (openssl, zlib, libcurl,
libuuid, libpulse). On a shared HPC system, some of these may only exist
as versioned environment modules (e.g. `cmake/3.26.3_gcc_9.5.0`, not bare
`cmake`) rather than on `PATH` by default, and if `svn` isn't available as
a module or system package at all, a non-root user can't add it —
`build_crux_from_source.sh` can't do anything about that case; ask your
HPC's support team. It can also take 30-90+ minutes with continuous
network access, so run it as its own multi-hour SLURM job — **not**
inside the short interactive session used for normal setup:
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

Tiers 2-4 all end up using crux **4.1/4.2** instead of **4.3**. That's
fine for this pipeline specifically: `05_percolator_to_moff.py` computes
retention time itself directly from the mzML files rather than relying on
tide-search's 4.3-only `retention time` column, so the moFF-bound output
is identical either way. (Tier 4's source build, if you ever need it, can
build any version including 4.3 — edit `CRUX_VERSION` in
`build_crux_from_source.sh` to match.)

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
