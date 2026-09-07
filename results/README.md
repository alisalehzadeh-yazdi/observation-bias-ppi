# Results

`scripts/run_single_transition.py` writes two CSVs per transition into this
folder (both ignored by git — see `.gitignore` — since they're regenerable):

- `jaccard_<database>_<organism>_v<from>-v<to>_n<runs>.csv` — one row per
  simulation run, one column per model, raw Jaccard similarity
  `J_i^(r)(m)` (cf. Fig. 2 / Fig. S1).
- `protein_intersections_<...>.csv` — same shape, each cell is the
  semicolon-joined list of proteins in the intersection of that run's
  simulated edges and the actually-added edges.

Bias scores `R_i(m)` (Figs. 3-5) are derived from the Jaccard CSVs via
`obsbias.metrics.bias_score` and are printed to stdout by the script; if you
want them persisted too, add a `--save-bias-scores` flag or compute them
from the checked-in Jaccard CSVs with a short follow-up script (not yet
included — see the top-level README's "Additional data needed" list).
