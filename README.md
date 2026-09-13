# Observation Bias in Evolving Protein Interaction Networks

Code accompanying:

> Salehzadeh-Yazdi A, Tuncbag N, Gursoy A, Keskin O, Hütt M-T. *A
> quantitative analysis of observation biases in protein interaction
> networks.* PLOS Computational Biology (submitted).

This repository quantifies whether successive releases of public
protein-protein interaction (PPI) databases (STRING, BioGRID, HIPPIE,
IntAct) expand by uniformly discovering unknown biology, or by
preferentially re-visiting already well-connected proteins. It does this by
comparing each empirically observed database update against six
edge-addition null/alternative models — RE, HCE, LCE, RNBE, RWE, CRWE (see
Table 1 of the manuscript) — using a Jaccard-overlap-based, cross-database
normalized bias score R_i(m).

## Repository layout

```
src/
  edge_addition_models.py   # the six models (Eq 5-11): weights + sampling
  bias_score.py              # Delta E_i, Jaccard, R_i(m) (Eq 1-3)
tests/
  test_models.py              # self-contained correctness checks, no pytest needed
scripts/
  example_pipeline.py         # end-to-end demo on synthetic data (no downloads needed)
data/
  README.md                   # what to download, and from where, for each database
results/
  README.md                   # what a real run writes here
docs/
  methods_to_code_map.md       # equation-by-equation index into the code
```

## Installation

```
pip install -r requirements.txt
```

Requires Python 3.9+.

## Quickstart

Run the test suite (no data download required, ~1 minute):

```
python3 tests/test_models.py
```

Run the end-to-end demo on a synthetic network (no data download required,
a few seconds):

```
python3 scripts/example_pipeline.py
```

Both scripts can be run from anywhere inside the repo; they locate `src/`
relative to their own file location.

## Using this on real data

`src/edge_addition_models.py` and `src/bias_score.py` operate on plain
`networkx.Graph` objects, so plugging in a real database transition is:

```python
import sys; sys.path.insert(0, "src")
import networkx as nx
from edge_addition_models import MODELS, IndexedGraph, precompute_node_scores, sample_new_edges
from bias_score import observed_delta_E, jaccard, bias_score

G_i  = nx.read_edgelist("data/processed/biogrid_2019_human.edgelist")
G_i1 = nx.read_edgelist("data/processed/biogrid_2022_human.edgelist")

delta_E = observed_delta_E(G_i, G_i1)
k = len(delta_E)
ig = IndexedGraph.build(G_i)

jaccards = {m: [] for m in MODELS}
for model in MODELS:
    scores = precompute_node_scores(ig, model) if model != "RNBE" else None
    for _ in range(100):  # n = 100 replicates, matching the manuscript
        run_scores = None if model == "RNBE" else scores
        sim = sample_new_edges(ig, k, model, node_scores=run_scores)
        jaccards[model].append(jaccard(sim, delta_E))

R = {m: bias_score(jaccards[m], jaccards["RE"]) for m in MODELS if m != "RE"}
print(R)
```

See `data/README.md` for where each database version comes from, and
`scripts/example_pipeline.py` for the fully worked version of this loop
(including how HCE/LCE/RWE/CRWE weights are precomputed once per
transition and reused across runs, while RNBE's node attribute is redrawn
every run, per the manuscript's Simulation Procedure).

## Correctness

Every equation this code implements (Eq 1–11 of the Methods) is checked in
`tests/test_models.py`: each model returns exactly the required number of
valid new edges; RE is verified not to be secretly degree-biased; HCE is
verified to systematically favor higher-degree endpoints than LCE; the
large-graph sampling scheme is checked against a closed-form theoretical
marginal (needed because the ~27,000-node BioGRID human network makes
enumerating the full candidate edge set impractical); and a full synthetic
Eq (1)–(3) pipeline is checked to recover a planted preferential-attachment
bias. Run `python3 tests/test_models.py` yourself — it prints a PASS line
for every check.

## Citation

See `CITATION.cff`. If you use this code, please cite the paper above.

## License

MIT — see `LICENSE`.
