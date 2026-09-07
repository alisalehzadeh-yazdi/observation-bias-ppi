# Observation Bias in Evolving Protein Interaction Networks

Code accompanying:

> Salehzadeh-Yazdi A, Tuncbag N, Gursoy A, Keskin O, Hütt M-T. *A
> quantitative analysis of observation biases in protein interaction
> networks.*

This repository quantifies whether successive releases of public
protein-protein interaction (PPI) databases (STRING, BioGRID, HIPPIE,
IntAct) expand by uniformly discovering unknown biology, or by
preferentially re-visiting already well-connected proteins. It does this by
comparing each empirically observed database update against edge-addition
null/alternative models, using a Jaccard-overlap-based, cross-database
normalized bias score.

**Before using this repo for anything citable**, read
`docs/code_consistency_review.md`. It documents exactly how this code was
checked against the manuscript's Methods section, three bugs that were
fixed during cleaning, and **two open questions that need a co-author
decision** (marked ‡ below) before the numbers this code produces can be
called the numbers behind Figs. 2-5.

## Repository layout

```
.
├── README.md                       <- you are here
├── LICENSE                         <- MIT (placeholder -- confirm with all co-authors)
├── CITATION.cff                    <- fill in venue/DOI once accepted
├── requirements.txt
├── .gitignore
├── docs/
│   └── code_consistency_review.md  <- code <-> manuscript diff, bugs fixed, open questions
├── src/obsbias/
│   ├── network_io.py                <- load + preprocess one database release
│   ├── growth_models.py             <- the six Table-1 models + 2 experimental ones
│   ├── simulate.py                  <- run_transition(): the core simulation loop
│   └── metrics.py                   <- Jaccard similarity + the bias score R_i(m)
├── scripts/
│   └── run_single_transition.py     <- CLI: reproduce one database/organism/version-pair
├── tests/
│   └── test_models.py               <- sanity tests (all passing; see below)
├── data/
│   └── README.md                    <- what raw files go where, download sources (‡ see below)
└── results/
    └── README.md                    <- output CSV schema
```

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .          # or just add src/ to PYTHONPATH, as scripts/ does

pytest tests/ -v          # 5/5 passing as of this cleaning pass

python scripts/run_single_transition.py \
    --early data/raw/string/Y4932.protein.physical.links.v11.5.txt \
    --late  data/raw/string/Y4932.protein.physical.links.v12.0.txt \
    --database STRING --organism yeast \
    --version-from 11.5 --version-to 12.0 \
    --threshold 400 --threshold-mode less \
    --n-runs 100 --outdir results/
```

This reproduces the specific run the original snippet computed (STRING
yeast, v11.5 → v12.0, `combined_score < 400`, one panel of Fig. S1),
using the paper's actual bias-score formula (see
`docs/code_consistency_review.md` §2.1 for why that matters) and without the
hardcoded/leftover names in the original script (§3).

## What changed from the pasted script (short version)

- Fixed the bias-score formula (`bias_score = mean(J(model)) / mean(J(RE))`,
  not `1 / mean(J(model))`).
- Fixed an edge-tuple-ordering bug that could leave a handful of
  already-existing edges in the "candidate" pool.
- Renamed all six models to the manuscript's own abbreviations
  (`RE, HCE, LCE, RNBE, RWE, CRWE`); the two models in the script that
  are **not** in the manuscript's Table 1 are kept, but clearly marked
  `*_experimental` and off by default (`--include-experimental` to turn
  them on).
- Removed hardcoded file names/dates and misleading variable names
  (`network_v9_1` that actually held v11.5 data); everything is now
  derived from CLI arguments.
- Fixed a stale print statement that echoed the wrong output filename.
- Turned a silent no-op (`if not enough candidates: print(...)`) into a
  raised `ValueError`.
- Replaced an O(|V|²) per-run rebuild in the "node-based random"
  experimental model with an O(|candidate edges|) one-time build.
- Added an optional `random_seed` for exact reproducibility.

Full details, including two items that need your (the authors') input
rather than a unilateral fix, are in `docs/code_consistency_review.md`.

## Additional data / materials needed before this repo is complete

Roughly in the order you'd want to tackle them:

1. **The real `process_file` / loaders for BioGRID, HIPPIE, and IntAct.**
   Only a STRING loader (`network_io.load_string_network`) is implemented
   here, reconstructed from the manuscript text because the original
   `process_file` function wasn't in the pasted snippet. Add
   `load_biogrid_network`, `load_hippie_network`, `load_intact_network` to
   `src/obsbias/network_io.py` following the same pattern (parse → drop
   self-loops/duplicates → optional LCC) once you have the original
   parsing code or the raw files' exact column layout.
2. **Confirm and resolve the two ‡ items in
   `docs/code_consistency_review.md`** (§2.2: are the two non-Table-1
   models intentional; §2.3: was the common-nodes-across-all-versions
   restriction actually used for the published figures, or only per-transition
   LCC). Whichever way these are resolved, update this repo's defaults
   (`restrict_to_common_nodes`, `include_experimental`) and the
   manuscript's Methods/Data sources paragraph to say so explicitly.
3. **The driver code for the full grid**, not just one transition: a
   script (e.g. `scripts/run_all_transitions.py`) that loops over all four
   databases × up to five organisms × all version-pairs × (for STRING) all
   four thresholding schemes, calling `run_transition` for each and
   writing one row per transition into a combined results table — this is
   what actually reproduces Figs. 2-5 and Figs. S1-S3, not the single-panel
   example script currently included.
4. **The figure-generation code** for Figs. 2-5 and S1-S3 (boxplots,
   the bias-score line/bar/heatmap plots). None of the plotting code from
   the original snippet was salvageable as-is (it plotted the wrong
   quantity — the mis-defined ratio, §2.1 — so it's not included here);
   add a `scripts/make_figures.py` that reads the combined results table
   from item 3 and reproduces each figure from the corrected `bias_score`.
5. **The edge-removal control analysis** (Fig. S3): the code for testing
   whether *removed* edges look random, referenced in the manuscript's
   Results ("In the curation process of STRING...") but not present in the
   pasted snippet at all.
6. **Raw data**, per `data/README.md`: download STRING/BioGRID/HIPPIE/IntAct
   releases, record exact source URLs/access dates in a `SOURCES.md` per
   database, and decide how large files are handled (excluded from git, as
   currently configured, vs. Git LFS vs. an external archive like Zenodo
   with a DOI linked from this README and the manuscript's Data
   Availability statement).
7. **Supplementary Table S1** (network statistics for all version
   transitions, referenced in the manuscript's Data sources paragraph) —
   either as a generated CSV (`results/table_s1_network_stats.csv`, one row
   per database/organism/version with node/edge counts before and after
   LCC) or committed as supplied by the authors.
8. **Environment pinning**: `requirements.txt` currently lists
   minimum versions only; once the full pipeline (item 3) has actually been
   run end-to-end, freeze exact versions (`pip freeze > requirements-lock.txt`)
   so the repo can reproduce the paper's exact numbers, not just numbers
   that are qualitatively consistent.
9. **A `setup.py` / `pyproject.toml`** so `pip install -e .` works as the
   README's Quick Start claims (not yet added — currently `scripts/`
   manually inserts `src/` onto `sys.path`, which works but is worth
   upgrading once the package layout above is finalized).
10. **Author contributions, funding, and a real Acknowledgments
    section** for the manuscript itself (currently "TBA") and a completed
    **Data and Code Availability statement** in the manuscript pointing at
    this repository's URL and, once minted, its Zenodo DOI (GitHub releases
    can be archived to Zenodo in a couple of clicks once the repo is public
    — see step-by-step below).

## Publishing this to GitHub

I don't have a connected GitHub account for this workspace, so I can't push
on your behalf, but everything above is ready to go from your machine:

```bash
cd <this folder>
git init
git add .
git commit -m "Initial cleaned release of observation-bias simulation code"
gh repo create <org-or-username>/observation-bias-ppi --public --source=. --remote=origin
git push -u origin main
```

(If you don't have the `gh` CLI, create the empty repo on github.com first,
then `git remote add origin <url>` and `git push -u origin main`.) Once
it's public, minting a Zenodo DOI for item 10 above is Settings → linking
the repo at zenodo.org/account/settings/github, then cutting a GitHub
Release.
