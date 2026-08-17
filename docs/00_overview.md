# 0. Overview: what this pipeline does and why

## The big picture

A mass spectrometer produces raw signal files (`.raw`) that record, for
thousands of tiny fragments of digested proteins ("peptides"), their mass
and abundance over time. Turning that raw signal into "here are the
proteins in my sample, and how much of each" takes four stages:

```
 .raw file                     the mass spectrometer's native output.
    |                          Unreadable by most downstream software.
    v
 [ThermoRawFileParser]         converts .raw -> the open standard .mzML
    |                          format that identification tools understand.
    v
 .mzML file
    |
    v
 [crux tide-search]            compares every observed spectrum against a
    |                          database of possible peptides (derived from
    |                          your protein FASTA) and proposes the best
    |                          matching peptide for each spectrum ("PSM" =
    |                          peptide-spectrum match). Includes many wrong
    |                          guesses at this stage — that's expected.
    v
 tide-search.txt (candidate matches, many false positives)
    |
    v
 [crux percolator]             uses machine learning + a "decoy" database
    |                          (reversed/shuffled fake proteins) to learn
    |                          which matches look real vs. random, and
    |                          keeps only matches below a chosen false
    |                          discovery rate (FDR), e.g. 10% or 1%.
    v
 percolator.target.psms.txt (confident, filtered peptide identifications)
    |
    v
 [moFF]                        goes back to the original .mzML signal and
    |                          measures the actual peak height/area (MS1
    |                          intensity) for each confidently identified
    |                          peptide in every sample — this is the
    |                          quantitative readout used to compare
    |                          protein abundance between samples/conditions.
    v
 quantified peptide/protein table  -> ready for statistics in R/Python/Excel
```

Everything here runs **locally on your own computer** — no cluster or
remote server needed. Just a Mac, Linux, or Windows machine with a few
GB of free disk space.

## Glossary

- **Peptide**: a short chunk of a protein produced by digesting the sample
  with an enzyme (usually trypsin) before it goes into the mass spec.
- **Spectrum / MS2 scan**: one fragmentation measurement of a single
  peptide, at one point in time (retention time).
- **PSM (peptide-spectrum match)**: a claim that "this spectrum = this
  peptide." Search engines like tide-search propose PSMs; most of the
  top-scoring PSMs from a single search are correct, but a meaningful
  fraction are wrong matches that happen to score well by chance.
- **Decoy database**: a fake protein database (e.g. every real sequence
  reversed) concatenated to the real ("target") database. Matches to decoy
  entries are known to be wrong by construction, so they tell you what
  "random noise" scores look like — this is what lets percolator estimate
  a false discovery rate.
- **FDR (false discovery rate)**: the estimated fraction of your *kept*
  matches that are actually wrong. "1% FDR" means about 1 in 100 of the
  matches you kept are expected to be false positives.
- **Retention time (RT)**: when, during the LC-MS run, a peptide eluted
  off the chromatography column and was measured. Needed by moFF to know
  where to look in the raw signal for a peptide's peak.
- **MS1 intensity / peak apex**: the height of a peptide's precursor-ion
  signal at its point of maximum abundance — this is what moFF extracts
  and reports as the quantitative value per peptide, per sample.
- **Label-free quantification (LFQ)**: comparing peptide/protein abundance
  across samples using raw signal intensity (as opposed to isotope-labeling
  methods like TMT/iTRAQ). moFF performs LFQ.
- **Match-between-runs (MBR)**: moFF can borrow a peptide's retention time
  from a run where it WAS confidently identified to look for (and quantify)
  the same peptide's signal in a run where it was missed by the ID search.
  This recovers more data points but should be used thoughtfully (see doc 6).

## The tools, in one line each

| Tool | Role | Docs |
|------|------|------|
| [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) | .raw → .mzML conversion | doc 3 |
| [crux](https://crux-toolkit.com) (`tide-search`) | peptide identification (search) | doc 4 |
| [crux](https://crux-toolkit.com) (`percolator`) | FDR filtering of identifications | doc 5 |
| [moFF](https://github.com/compomics/moFF) | label-free MS1 quantification | doc 6 |

## Why is there a "proteomics" conda environment if crux isn't in it?

crux and ThermoRawFileParser are downloaded directly as ready-to-run
binaries (doc 1) rather than installed via conda — they're not available
through conda on all platforms this course supports (Mac, Linux,
Windows), so a direct download keeps the install method identical
everywhere. The **`proteomics`** conda environment just holds Python +
`lxml` for the small script that bridges crux's output into moFF's input
format (doc 6). moFF gets its own separate **`moff`** environment because
it's unmaintained software pinned to Python 3.6 with old, exact versions
of numpy/pandas/scikit-learn — incompatible constraints that don't belong
in the same environment as everything else. Every script in this repo
tells you which environment (if any) it needs.

## What about Windows?

Every tool here (crux, ThermoRawFileParser, moFF) has a genuine, official
Windows build, so Windows users run this pipeline with the same steps
and same tool versions as Mac/Linux users — no virtual machine or WSL
required. The scripts live in a separate `windows/` folder (mirroring
`mac_linux/` one-for-one) so nothing depends on runtime OS-detection, and
they run inside **Git Bash** (from "Git for Windows"), since the scripts
are bash. See [doc 1b](01b_windows_setup.md).

Next: [1. Installing the software](01_install_software.md)
