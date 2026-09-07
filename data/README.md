# Data

Raw database releases are **not** committed to this repository (see
`.gitignore`) — they are large, versioned externally, and in some cases
(BioGRID, IntAct) redistribution has its own license terms. This folder
documents exactly what each script expects and where to get it.

## Expected layout

```
data/
  raw/
    string/    Y4932.protein.physical.links.v9.1.txt ... v12.0.txt
               9606.protein.physical.links.vX.Y.txt  (human)
               <taxon>.protein.physical.links.vX.Y.txt (E. coli, Drosophila, mouse)
    biogrid/   BIOGRID-ORGANISM-<Organism>-<version>.tab3.txt   (one per organism/version)
    hippie/    hippie_v2.X.txt                                  (human only)
    intact/    intact_<year>.txt (or the PSI-MI TAB export for the relevant release)
  processed/
    <database>_<organism>_v<version>.graphml   (optional cache of loaded/cleaned graphs)
```

## Download sources

| Database | Where | Versions used in the paper |
|---|---|---|
| STRING | https://string-db.org/cgi/download | v9.1, v10, v10.5, v11, v11.5, v12 |
| BioGRID | https://downloads.thebiogrid.org/BioGRID | 2013, 2016, 2019, 2022, 2025 releases, for *E. coli*, *S. cerevisiae*, *D. melanogaster*, *M. musculus*, *H. sapiens* |
| HIPPIE | http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/ | 2010, 2014, 2020, 2026 (human only) |
| IntAct | https://www.ebi.ac.uk/intact/ | 2011, 2013, 2016, 2019, 2022 (human only) |

## Step-by-step to add real data

1. Download each release listed above into `data/raw/<database>/` using
   each database's own version-archive download page (most keep old
   releases at a stable URL; STRING and BioGRID both version their FTP/HTTP
   download paths).
2. Record the exact download URL and access date for each file in
   `data/raw/<database>/SOURCES.md` (create one per database folder) — this
   is what a reader will need to verify or refresh the dataset later, and
   what a "Data Availability" statement in the manuscript should point to.
3. If a file's raw format differs from STRING's `protein1 protein2
   combined_score` layout (BioGRID, HIPPIE, and IntAct all do), add a
   loader to `src/obsbias/network_io.py` following the pattern of
   `load_string_network` — see the top-level README's "Additional data
   needed" section, item 2.
4. Do not commit raw files larger than GitHub's soft 50MB / hard 100MB
   per-file limits; for large releases (BioGRID/IntAct full dumps can
   exceed this) either keep them out of git entirely (as configured here)
   or use Git LFS / an external archive (e.g., Zenodo) and link to it from
   this file.
