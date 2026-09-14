# Methods equations → code

A direct index from each numbered equation in the manuscript's Methods
section to the function that implements it, so a reviewer (or future you)
can check the code against the paper line by line without hunting for it.

| Equation | What it defines | Code |
|---|---|---|
| Eq (1): ΔE_i = E_{i+1} \\ E_i | The empirical update between two consecutive database versions, restricted to node pairs present in both versions | `src/bias_score.py` → `observed_delta_E(G_i, G_i1)` |
| Eq (2): J_i^(r)(m) | Jaccard similarity between one simulated update and the real update | `src/bias_score.py` → `jaccard(A, B)` |
| Eq (3): R_i(m) | Bias score: mean Jaccard for model *m* normalized by mean Jaccard for RE | `src/bias_score.py` → `bias_score(J_model_runs, J_RE_runs)` |
| Eq (4): p_m(u,v) | Sampling probability of a candidate edge, proportional to its model weight | `src/edge_addition_models.py` → `sample_new_edges` (both the small-graph exact path and the large-graph scalable path implement this same target distribution; see the module docstring for why two code paths exist and how they're validated against each other) |
| Eq (5): RE | Uniform sampling weight | `precompute_node_scores(..., "RE")` returns `None` (no per-node score needed); `_generate_candidate_pairs_uniform` |
| Eq (6): HCE, w = k_u + k_v | Degree-sum weight | `precompute_node_scores(..., "HCE")` |
| Eq (7): LCE, w = k_max/(k_u+1) + k_max/(k_v+1) | Inverse-degree-sum weight | `precompute_node_scores(..., "LCE")` |
| Eq (8): RNBE, w = r_u · r_v | Product of per-node Uniform(0,1) attributes, redrawn every simulation run | `precompute_node_scores(..., "RNBE")`; see `sample_new_edges`'s handling of `node_scores=None` for the "redraw every run" behavior |
| Eq (9): RWE, w = π_u + π_v | PageRank-sum weight (α = 0.85, networkx implementation) | `precompute_node_scores(..., "RWE")` |
| Eq (10): s_v | Compact random walk score (N_walks=5000, L=100, p_r=0.5) | `compact_random_walk_scores(...)` |
| Eq (11): CRWE, w = s_u + s_v | Compact-random-walk-score-sum weight | `precompute_node_scores(..., "CRWE")` |

## Simulation Procedure → code

The manuscript's "Simulation Procedure" paragraph specifies that RE, HCE,
LCE, RWE, and CRWE weights are precomputed once per version transition and
reused across all `n = 100` replicates, while RNBE's node attribute is
redrawn independently at the start of every replicate. This is implemented
by calling `precompute_node_scores` once per transition for the first five
models (its return value is then passed into `sample_new_edges` on every
replicate), and by passing `node_scores=None` into `sample_new_edges` for
RNBE on every replicate, which makes it draw a fresh attribute each time.
`scripts/example_pipeline.py` shows this pattern end to end.
