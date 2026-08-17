# 3. Converting .raw files to .mzML

Thermo mass spectrometers save data in a proprietary `.raw` format that
identification software can't read directly. The first step is converting
each `.raw` file to the open **mzML** format using ThermoRawFileParser.

## 3.1 Organize your raw files

Put all the `.raw` files for a project in one folder on your computer,
e.g.:

```bash
mkdir -p ~/my_project/raw
# copy or move your .raw files into ~/my_project/raw/
```

## 3.2 Run the conversion

```bash
cd proteomics_pipeline_course/mac_linux/scripts   # or windows/scripts on Windows
./01_convert_raw_to_mzml.sh ~/my_project/raw ~/my_project/mzML
```

This calls ThermoRawFileParser (installed by doc 1 / 1b) once per `.raw`
file and writes one matching `.mzML` file per input, e.g. `1.raw` →
`1.mzML`. Keep the same base filenames throughout the pipeline — moFF
(doc 6) later depends on filenames matching up between steps.

For many files or large files, this step can take a while (it's I/O and
CPU bound) — let it run in the background or overnight if needed; there's
no need to babysit it.

## 3.3 Sanity-check the output

```bash
ls -lh ~/my_project/mzML
```

You should see one `.mzML` file per `.raw` file, each with a substantial
size (tens to low hundreds of MB is typical — an mzML is usually *larger*
than the source .raw file, since it stores data as text/XML).

## 3.4 If a file fails with "still being acquired"

```
ERROR RAW file cannot be processed since it is still being acquired
```

This means that specific `.raw` file's internal "acquisition complete" flag
was never set — usually because the run was stopped early / aborted by the
instrument operator, or the file was copied off the instrument before
Xcalibur fully closed it. It is a property of that particular file, not a
bug in the conversion step (other files converting fine confirms this). See
[doc 7](07_troubleshooting.md) for how to confirm and what your options are.

Next: [4. Peptide identification with tide-search](04_tide_search.md)
