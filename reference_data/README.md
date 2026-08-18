# reference_data/

Shared reference files used automatically by both tracks' database-prep
step (`02_build_database_index.sh`) — you don't need to do anything with
these files yourself.

## `common_contaminants.fasta`

The classic cRAP ("common Repository of Adventitious Proteins") database
from the Global Proteome Machine organisation — 116 sequences: human
keratins, trypsin, BSA, and other proteins that routinely turn up in mass
spec samples regardless of what you're actually studying, from lab
handling and sample prep. `02_build_database_index.sh` appends this to
your database automatically before indexing (see
[docs/02_prepare_database.md](../docs/02_prepare_database.md)) — without
it, contaminant peptides get misassigned to whatever real protein in your
database happens to look similar, which is a real source of incorrect
identifications, not just a cosmetic issue.

Source: mirrored from
[Zenodo record 10.5281/zenodo.15115102](https://zenodo.org/records/15115102)
(`crap_gpm.fasta`), licensed CC-BY-4.0. Original: <https://www.thegpm.org/crap/>.
