# Data

Raw database releases are **not** committed to this repository (see
`.gitignore`) — they are large, versioned externally, and in some cases
redistribution has its own license terms. This folder documents exactly
what each script expects and where to get it.

## Download sources

| Database | Where | Versions used in the paper |
|---|---|---|
| STRING | https://string-db.org/cgi/download | v9.1, v10, v10.5, v11, v11.5, v12 |
| BioGRID | https://downloads.thebiogrid.org/BioGRID | 2013, 2016, 2019, 2022, 2025 releases, for *E. coli*, *S. cerevisiae*, *D. melanogaster*, *M. musculus*, *H. sapiens* |
| HIPPIE | http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/ | 2011, 2013, 2016, 2019, 2022 (human only) |
| IntAct | https://www.ebi.ac.uk/intact/ | 2010, 2014, 2020, 2026 (human only) |

Exact node/edge counts for every version and organism above are reported in
Supplementary Table S1 of the manuscript — use those numbers to sanity-check
your own downloads before running the analysis.

## Step-by-step to add real data

1. Download each release listed above into `data/raw/<database>/`, using
   each database's own version-archive download page (STRING and BioGRID
   both version their FTP/HTTP download paths; HIPPIE and IntAct keep
   dated release archives).
2. Record the exact download URL and access date for each file in
   `data/raw/<database>/SOURCES.md` (create one per database folder) —
   this is what lets someone else reproduce the exact input you used,
   since these databases update in place and don't always keep old
   versions at a stable, permanently-dated URL.
3. Preprocess each raw file into a plain undirected edge list before
   handing it to `src/edge_addition_models.py`: remove self-loops and
   duplicate edges, and restrict to the largest connected component, per
   the manuscript's Data sources paragraph. For STRING, additionally
   filter by `combined_score` at the threshold you're reproducing (all
   interactions, `combined_score > 400`, `combined_score > 900`, or
   `combined_score < 400`).
4. Save the preprocessed network for each version under
   `data/processed/<database>_<version>_<organism>.edgelist` (or your own
   naming scheme) — `networkx.write_edgelist` / `networkx.read_edgelist`
   round-trip this format directly, matching the loading snippet in the
   top-level README.

## Why raw dumps aren't committed here

BioGRID, STRING, HIPPIE, and IntAct each publish their own reuse terms, and
the raw files for the version range used in this study (five databases ×
up to six versions × up to five organisms) are well beyond a size that
belongs in a git repository. Pointing to the official, versioned source is
also what keeps this reproducible if any of these databases changes its
raw-file format in the future — the preprocessing step above is what
absorbs that, rather than a stale binary blob.
