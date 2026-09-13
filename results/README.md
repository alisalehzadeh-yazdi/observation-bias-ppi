# Results

This folder is where output from running the analysis on real data belongs
— e.g. per-transition CSVs of Jaccard similarities (Eq 2) and bias scores
(Eq 3), one file per database/organism/version-transition, the raw
material behind Figs. 2-5 and S1-S3.

Nothing is committed here yet beyond this README: outputs are regenerable
from `data/` plus the code in `src/`, so keeping them out of git (see
`.gitignore`) avoids committing large, derived, easily-regenerated files.
As the analysis scripts for real data are added, this file should be
updated to describe their exact output format (columns, naming
convention), the way `data/README.md` documents the input side.
