# 8. Running this pipeline on Mizzou Hellbender (optional)

Everything in docs 0–7 runs entirely on your own computer — you don't need
this doc at all for that. This doc is for **larger, real studies** where
you want more CPU/memory than your laptop has, using the University of
Missouri System's **Hellbender** HPC cluster.

> Cluster details (hostnames, quotas, partition names) change over time.
> This doc reflects Hellbender's setup as documented at the time of
> writing — if anything here doesn't match what you see, trust the
> official docs at <https://docs.itrss.umsystem.edu/pub/hpc/hellbender>
> over this file, and let your instructor/PI know so it can be updated.

## 8.1 Get an account

Hellbender requires a **faculty sponsor** (your PI is listed as the
sponsor/PI on the request). Submit the request form here:
<https://tdx.umsystem.edu/TDClient/36/DoIT/Requests/ServiceOfferingDet?ID=1041>

For general questions, contact `itrss-support@umsystem.edu`.

## 8.2 Log in

**Mac/Linux/Windows (Git Bash)** — same command everywhere:

```bash
ssh yourusername@hellbender-login.rnet.missouri.edu
```

Use your campus SSO username and password. **Windows users can use Git
Bash** (see doc 1b) instead of a separate SSH client — `ssh` just works
the same as Mac/Linux. PuTTY is the traditional recommendation if you'd
rather use that.

**Off-campus, without VPN:** password login is disabled off-campus by
default — you need public-key authentication instead:
```bash
ssh-keygen -t ed25519 -C "you@missouri.edu"
```
Then email your **public** key (the `.pub` file, never the private one)
to `itrss-support@umsystem.edu` with the subject line **"Hellbender
Public Key Access"**. Easiest alternative: connect to the campus VPN
first, then plain password SSH works as above.

## 8.3 ⚠️ Never run this pipeline on the login node

The node you land on after `ssh` is a shared **login node** — for file
management and submitting jobs only, never for actual computation. All
the real pipeline steps (conversion, search, percolator, moFF) must run
as **SLURM batch jobs** (below) or in an interactive session:

```bash
salloc --partition=interactive --time=1:00:00
```

## 8.4 Storage layout

| Location | Quota | Use for |
|---|---|---|
| `$HOME` (`/home/$USER`) | 50 GB | Your code/scripts (this repo) |
| `$HOME/data` | 500 GB | Raw data, mzML files, pipeline output (`work/`) |
| `/local/scratch` | 1.6–3.2 TB, per-node, temporary | Not used by this pipeline |

**Nothing here is backed up.** Copy results you want to keep somewhere
else when you're done.

Put the actual repository (and your raw data) under `$HOME/data/`, not
directly in `$HOME` — the scripts below assume
`$HOME/data/proteomics_pipeline_course`.

```bash
mkdir -p $HOME/data
cd $HOME/data
# from your local machine, in a separate terminal:
#   rsync -av /local/path/to/proteomics_pipeline_course/ yourusername@hellbender-login.rnet.missouri.edu:~/data/proteomics_pipeline_course/
```

## 8.5 conda on Hellbender

Don't run `00_install_conda.sh` here — conda is already available as a
module:

```bash
module load miniconda3
```

By default conda stores environments/packages under `$HOME/.conda`,
which can quietly eat into your 50 GB `$HOME` quota (crux/moFF's
dependencies add up). Point it at your larger `$HOME/data` allocation
instead, **once**, before creating any environments:

```bash
mkdir -p $HOME/data/miniconda/envs $HOME/data/miniconda/pkgs
cat >> ~/.condarc << 'EOF'
envs_dirs:
  - /home/YOURUSERNAME/data/miniconda/envs
pkgs_dirs:
  - /home/YOURUSERNAME/data/miniconda/pkgs
EOF
```

(Replace `YOURUSERNAME` with your actual username — `.condarc` doesn't
expand `$HOME` or `$USER`.)

Then run the normal setup (in an interactive session, doc 8.3 — not the
login node):

```bash
salloc --partition=interactive --time=1:00:00
cd $HOME/data/proteomics_pipeline_course/mac_linux/scripts
module load miniconda3
./00_setup_environment.sh
exit    # leave the interactive session when done
```

### Why crux needs a fallback on Hellbender specifically

`00_setup_environment.sh` normally installs crux by downloading a
prebuilt binary from GitHub — fast, and no conda solving needed. That
binary needs glibc 2.29+, but Hellbender's OS (AlmaLinux 8.9) ships glibc
2.28 (check your own node with `ldd --version`), one minor version short.
This is confirmed working end-to-end on Hellbender as of this writing,
not just a theoretical fix:

1. The script checks your system's glibc version *before* downloading
   anything (via `ldd --version`), so on Hellbender it skips the doomed
   direct download entirely rather than wasting time on it.
2. It installs crux from conda instead, pinned to a specific bioconda
   build (`crux-toolkit=4.1=h503566f_3`) verified — by inspecting the
   actual GLIBC/GLIBCXX symbol versions compiled into that binary — to
   need only glibc 2.14, comfortably under Hellbender's 2.28. (Bioconda's
   *newest* crux build, 4.2, has the same glibc-2.29 problem as the direct
   download and does **not** work here; 4.1 does. Both were released the
   same day with an identical CLI, so this doesn't cost you anything
   functionally.)
3. It replaces `software/crux/bin/crux` with a small wrapper script that
   activates that conda environment and runs the real `crux` — and that
   wrapper loads the `miniconda3` module itself if `conda` isn't already
   on `PATH` in whatever shell calls it. So **you never need to manually
   run `module load miniconda3` / `conda activate crux-fallback` before
   using crux** — just call `software/crux/bin/crux ...` (or run the
   pipeline scripts, which all call it the same way) directly, in any
   shell, including inside a SLURM job. `00_setup_environment.sh`'s final
   "Setup complete" message tells you which backend ended up installed.

If setup ever gets interrupted partway, or you want to force a clean
reinstall of just crux, delete it and re-run:
```bash
rm -f $HOME/data/proteomics_pipeline_course/mac_linux/software/crux/bin/crux
./00_setup_environment.sh
```

On the rare chance the conda fallback fails for some other reason on your
account, the script automatically tries running crux inside a
Singularity/Apptainer container next (needs Singularity/Apptainer enabled
for your account — Hellbender describes this as "limited support ... for
advanced users"; if missing, email `itrss-support@umsystem.edu`), and as
a last resort prints instructions for compiling crux from source (a real,
slow build, and on Hellbender specifically may be blocked entirely if
`svn` isn't available even as a module — see the SLURM template in the
next section). See
[doc 7](07_troubleshooting.md#crux-libm-so6-version-glibc_229-not-found-common-on-hellbender--older-hpc-linux)
for the full details on all of this.

## 8.6 Using the SLURM templates

`mac_linux/scripts/slurm/` has one `.slurm` file per pipeline step, plus
`submit_all.sh` which submits all of them as a **dependency chain** (each
step only starts after the previous one succeeds). There's also a
standalone `build_crux_from_source.slurm`, not part of that chain — only
needed if doc 8.5's crux fallbacks both failed (see above).

```bash
cd $HOME/data/proteomics_pipeline_course/mac_linux/scripts/slurm
mkdir -p logs

# 1. First, edit the path variables near the top of EACH .slurm file to
#    match your project (they default to $HOME/data/proteomics_pipeline_course
#    and the bundled example data).

# 2. Then submit the whole chain:
./submit_all.sh
```

Watch progress with `squeue --me`. Check `logs/*.out` for each step's
output — that's where errors show up if a job fails.

You can also submit steps one at a time while learning the pipeline
(recommended the first time, so you can check each output before moving
on):
```bash
sbatch 01_convert.slurm
squeue --me                  # wait for it to finish
sbatch 02_index.slurm
# ...and so on
```

### If `sbatch` rejects your job asking for an account

Each `.slurm` file has a commented-out `#SBATCH --account=your_pi_account`
line. If submission fails complaining about a missing/invalid account,
ask your PI what Hellbender account name to use, uncomment that line, and
fill it in.

### Choosing resources

The templates have reasonable starting defaults for the small example
dataset. For a real study:

- `01_convert` (RAW→mzML): scales with number and size of raw files.
- `02_index`: fast regardless of size, rarely needs adjustment.
- `03_search`: the slowest step — scales with database size × number of
  spectra. A few thousand proteins × several runs can take **many hours**;
  give it generous `--time`. The `general` partition caps individual jobs
  at 2 days — if you need longer, ask your PI whether your group has
  access to a priority/investor partition (up to 28 days).
- `04_percolator`: usually fast (minutes).
- `05_moff`: scales with number of confidently identified peptides ×
  number of runs; can also take a while for larger studies.

If a job runs out of time or memory, check the log file for the specific
error, increase the relevant `#SBATCH` value, and resubmit — nothing about
re-running these scripts is destructive as long as you don't delete your
input data.
