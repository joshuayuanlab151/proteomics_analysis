# 2. Preparing your protein database

Every peptide identification search needs a **protein FASTA database**: a
text file listing every protein sequence the search engine is allowed to
consider as a possible match. Only peptides that come from proteins in this
file can be identified — so the database has to be relevant to your sample
(e.g. the correct species' proteome).

## 2.1 FASTA format

```
>sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGH
GKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEF
TPAVHASLDKFLASVSTVLTSKYR
```

Each entry starts with `>` followed by an identifier/description line, then
the amino acid sequence (usually wrapped to 60-80 characters per line).

Count how many proteins are in a FASTA file:

```bash
grep -c '>' your_database.fasta
```

## 2.2 Where to get a database

- **Human / common model organisms**: [UniProt](https://www.uniprot.org/)
  (download the reviewed "Swiss-Prot" proteome for your organism).
- **Plant genomes**: [Phytozome](https://phytozome-next.jgi.doe.gov/)
- **Fungal genomes**: [MycoCosm](https://mycocosm.jgi.doe.gov/mycocosm/home)
- Ask your PI/instructor if there's a lab-standard database for your project
  (e.g. a custom genome annotation).

## 2.3 (Optional) remove redundant sequences

If your database has many near-duplicate sequences (common with some genome
annotations), you can collapse them with CD-HIT before searching, which
speeds up the search and simplifies protein-level results:

```bash
conda install -c bioconda cd-hit
cd-hit -c 1.0 -n 5 -i your_input_proteome.fasta -o your_proteome_dedup.fasta -d 0
```

`-c 1.0` keeps only sequences that are exactly identical to each other
merged; lower it (e.g. `-c 0.95`) to also collapse near-identical sequences.

## 2.4 (Optional) add common contaminants

Lab contaminants (human keratin, trypsin itself, etc.) commonly appear in
samples and can otherwise get misassigned to real proteins. Many labs
append a standard contaminant FASTA (e.g. the
[MaxQuant contaminants list](https://www.maxquant.org/) or the
[Global Proteome Machine common contaminants](https://www.thegpm.org/crap/))
to their database:

```bash
cat your_proteome.fasta contaminants.fasta > combined_database.fasta
```

## 2.5 Decoys — handled automatically, don't do this yourself

Older versions of this pipeline had you build a separate target-decoy FASTA
by hand (e.g. with OpenMS `DecoyDatabase`). **You do not need to do this
anymore** — `crux tide-index` (the next step, doc 4) generates decoy
sequences internally and builds the target+decoy search index in one
command. Just point it at your real (target-only) FASTA.

## 2.6 The example database

`example_data/PSM_BL.fasta` is a small (22-protein) test database included
with this repo so you can run the whole pipeline end-to-end quickly before
using it on real, larger data. It is **not** a real study database — don't
use it for actual research.

Next: [3. Converting RAW files to mzML](03_convert_raw_to_mzml.md)
