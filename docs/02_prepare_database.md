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

## 2.3 Remove redundant sequences (required)

Exact-duplicate sequences in a database aren't harmless — they make the
search engine and downstream tools double-count/misreport peptides that
are shared across the duplicate entries, which distorts protein-level
results. `02_build_database_index.sh` handles this **automatically**, for
every run, before indexing: it collapses any sequences that are exactly
identical to each other down to one entry each, using a small built-in
`awk` step (no separate tool to install, and it works the same on the
Windows track). You don't need to do anything yourself — the script logs
how many duplicates it removed each time it runs.

If your database has many *near*-identical (not exactly identical)
sequences — common with some genome annotations, less common with
UniProt/Swiss-Prot downloads — the built-in step won't catch those, since
that requires real sequence-similarity clustering. On Mac/Linux only, you
can optionally run [CD-HIT](https://github.com/weizhongli/cdhit) yourself
first for that case (CD-HIT has no Windows build):

```bash
conda install -c bioconda cd-hit
cd-hit -c 0.95 -n 5 -i your_input_proteome.fasta -o your_proteome_dedup.fasta -d 0
```

then pass `your_proteome_dedup.fasta` to `02_build_database_index.sh`
instead of your original file — it'll still do the exact-duplicate pass
and contaminant append on top of that.

## 2.4 Add common contaminants (required)

Lab contaminants (human keratins, trypsin itself, BSA, etc.) are present
in essentially every real mass spec sample, regardless of what you're
actually studying, from sample handling and prep. If they're not in your
search database, those peptides don't just go unidentified — they get
**misassigned to whatever real protein in your database happens to look
similar**, which is a real source of incorrect identifications, not a
cosmetic gap.

`02_build_database_index.sh` handles this **automatically** too: it
appends the bundled
[common proteomics comtaminations](https://zenodo.org/records/15115102)
(the standard 116-sequence cRAP database) to your database before
indexing, every time. Again, nothing for you to do — the script logs how
many contaminant sequences it added.

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

Like any other database, it still goes through 2.3/2.4 automatically —
`02_build_database_index.sh` actually indexes 22 + 116 = 138 sequences for
it (the 22 real proteins plus the full contaminants list), even though the
file on disk has 22. If you're ever counting proteins with
`grep -c '>' your_database.fasta` (2.1) to sanity-check something, count
against the `..._prepared.fasta` file next to your index output, not your
original input file, to see what was actually searched.

Next: [3. Converting RAW files to mzML](03_convert_raw_to_mzml.md)
