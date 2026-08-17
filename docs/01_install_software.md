# 1. Installing the software (Mac/Linux — do this once)

**Using Windows?** Skip this doc and read [doc 1b](01b_windows_setup.md)
instead — it covers the equivalent Windows setup.

Everything below assumes you have this whole `proteomics_pipeline_course`
folder on your computer already (unzip it or `git clone` it wherever you
like). All commands below use the `mac_linux/` track's scripts.

## 1.1 Do you already have conda?

Conda (via Miniconda or Anaconda) is what we use for moFF's Python
environment and the small helper script that bridges crux's output to
moFF's input format.

```bash
conda --version
```

- If this prints a version number, skip to 1.2.
- If you get "command not found", install it with the provided script:
  ```bash
  cd proteomics_pipeline_course/mac_linux/scripts
  ./00_install_conda.sh
  ```
  This installs Miniconda into `$HOME/miniconda3` and wires it up for
  bash. **Close and reopen your terminal** afterwards, then confirm
  `conda --version` works before continuing.

## 1.2 Run the setup script

```bash
cd proteomics_pipeline_course/mac_linux/scripts
./00_setup_environment.sh
```

This downloads two tools directly (official builds, not through conda —
crux and ThermoRawFileParser aren't available via conda on every platform
this course supports) and creates two conda environments:

1. **crux** — downloaded to `mac_linux/software/crux/` (bundles both
   `tide-search` and `percolator` — you do not need to install percolator
   separately).
2. **ThermoRawFileParser** — downloaded to
   `mac_linux/software/ThermoRawFileParser/` (self-contained; no separate
   Mono or .NET install needed).
3. **`proteomics`** conda environment — just Python + the `lxml` library,
   used by the small script that converts crux's output into moFF's input
   format (doc 6).
4. **`moff`** conda environment — the exact old package versions moFF
   needs, plus downloads and patches the moFF source code itself into
   `mac_linux/software/moFF/` (see doc 7 for what the patch fixes and why).

This takes several minutes (mostly downloading). You only need to do this
once per computer.

### ⚠️ Do NOT `pip install moFF`

PyPI has a package literally named `moFF` (and `MOFF`) — but it is a
**completely different tool** for predicting CRISPR/Cas9 off-target effects,
from an unrelated group at MD Anderson. It is not related to mass
spectrometry at all. If you `pip install moFF` you will get the wrong
software with no obvious error. Always install moFF via the setup script
above (which clones the real tool from
[compomics/moFF](https://github.com/compomics/moFF) on GitHub).

## 1.3 Verify the install

From inside `mac_linux/scripts/`:

```bash
../software/crux/bin/crux version
../software/ThermoRawFileParser/ThermoRawFileParser --version

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate moff
python ../software/moFF/moff_all.py -h
```

If all three commands print version/help text without errors, you're set.
See [doc 7](07_troubleshooting.md) if something fails — a couple of known,
easy-to-fix issues (e.g. macOS refusing to run downloaded binaries) are
documented there.

Next: [2. Preparing the protein database](02_prepare_database.md)
