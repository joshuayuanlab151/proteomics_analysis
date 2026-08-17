# 1b. Windows-specific setup

Everything in this pipeline (crux, ThermoRawFileParser, moFF) has genuine,
official Windows builds — you do **not** need a Linux virtual machine or
WSL to run this on Windows. What you do need is a **bash shell**, because
all the scripts in `windows/scripts/` are bash (`.sh`), not `.bat`/
PowerShell. The easiest way to get one on Windows is **Git Bash**, which
also happens to install `git` (needed anyway, to fetch moFF's source).

> **This Windows path has not been tested on a real Windows machine** while
> building this pipeline (only Mac, and implicitly Linux). The underlying
> tools (crux, ThermoRawFileParser) ship official Windows binaries and
> `windows/scripts/` downloads and invokes them directly, so it is
> expected to work — but if you hit something that doesn't, please report
> exactly what command and error you saw. If the automated steps below
> ever fail partway, each one has a documented manual fallback — look for
> the ⚠️ notes.

## 1. Install Git for Windows (gives you Git Bash)

Download and run the installer from <https://git-scm.com/download/win>.
Default options are fine. This adds a **"Git Bash"** entry to your Start
Menu — that's the shell you'll use for everything in this pipeline from
now on (never Command Prompt or PowerShell).

## 2. Get this repository onto your machine

If you were given a zip file, extract it anywhere convenient (e.g.
`C:\Users\you\proteomics_pipeline_course`) and open **Git Bash**, then
`cd` into it:

```bash
cd /c/Users/you/proteomics_pipeline_course
```

(Git Bash uses `/c/...` instead of `C:\...` — this is normal. If you ever
need to type a Windows path by hand, e.g. copied from File Explorer,
convert `C:\Users\you\...` to `/c/Users/you/...` first.)

## 3. Install conda

```bash
conda --version
```

If that fails with "command not found":

```bash
cd windows/scripts
./00_install_conda.sh
```

This downloads and silently installs Miniconda, then wires it up for Git
Bash. **Close and reopen Git Bash** afterwards, then confirm
`conda --version` works.

> ⚠️ **If `00_install_conda.sh` fails**, or warns about a space in your
> Windows username/path, install Miniconda manually instead:
> 1. Download <https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe>
> 2. Double-click it and click through the installer normally (default
>    options are fine; "Just Me" install is recommended).
> 3. Re-run `./00_install_conda.sh` anyway afterward — it will detect the
>    existing install and just do the Git Bash wiring step.

## 4. Run the pipeline software setup script

```bash
cd windows/scripts    # if not already there
./00_setup_environment.sh
```

This downloads the Windows builds of crux and ThermoRawFileParser and sets
up the same two conda environments (`proteomics`, `moff`) as the Mac/Linux
track.

### Expect a SmartScreen warning

The first time `ThermoRawFileParser.exe` or `crux.exe` actually runs,
Windows may show a **"Windows protected your PC"** SmartScreen prompt,
because these are downloaded executables without a Microsoft-recognized
signature. This is expected — click **"More info"** then **"Run anyway"**.
You should only need to do this once per executable.

## 5. Verify the install

From inside `windows/scripts/`:

```bash
../software/crux/bin/crux.exe version
../software/ThermoRawFileParser/ThermoRawFileParser.exe --version

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate moff
python ../software/moFF/moff_all.py -h
```

If all three commands print version/help text without errors, you're set.

## 6. Everything else

From here on, follow [docs/02](02_prepare_database.md) through
[docs/07](07_troubleshooting.md) — wherever they say `scripts/whatever.sh`,
use `windows/scripts/whatever.sh` (the script names and arguments are
identical to Mac/Linux, only the location differs).

Next: [2. Preparing the protein database](02_prepare_database.md)
