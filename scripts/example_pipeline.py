"""
example_pipeline.py

A minimal, runnable, end-to-end demonstration of the full analysis
pipeline (Eq 1-11 of Methods) on synthetic data, so a reader can see how
the pieces fit together without needing any of the real BioGRID / STRING /
HIPPIE / IntAct downloads. Swap `G_i` / `G_i1` below for two real
consecutive-version networkx graphs (see download scripts) to run this on
actual data -- everything past that point is identical.

Run with:  python3 example_pipeline.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import networkx as nx

from edge_addition_models import (
    MODELS,
    IndexedGraph,
    precompute_node_scores,
    sample_new_edges,
)
from bias_score import observed_delta_E, jaccard, bias_score

N_SIMULATION_RUNS = 100  # matches "n = 100 independent simulation runs" in Methods


def load_toy_transition(seed=0):
    """Stand-in for loading two consecutive real database versions.
    Replace this with, e.g.:

        G_i  = load_network("data/processed/biogrid_2019_human.edgelist")
        G_i1 = load_network("data/processed/biogrid_2022_human.edgelist")

    using whatever edge-list format your `src/download/` scripts produce.
    """
    rng = np.random.default_rng(seed)
    G_i = nx.barabasi_albert_graph(300, 3, seed=seed)

    # Build a synthetic "next version" by adding edges with a genuine
    # (but not perfect) preference for high-degree endpoints, so the demo
    # has a real, recoverable bias signature rather than pure noise.
    G_i1 = G_i.copy()
    deg = dict(G_i.degree())
    nodes = list(G_i.nodes())
    weights = np.array([deg[v] + 1.0 for v in nodes])
    weights /= weights.sum()
    target_new_edges = 120
    added = 0
    while added < target_new_edges:
        u, v = rng.choice(nodes, size=2, replace=False, p=weights)
        if not G_i1.has_edge(u, v):
            G_i1.add_edge(u, v)
            added += 1
    return G_i, G_i1


def run_transition(G_i: nx.Graph, G_i1: nx.Graph, n_runs: int = N_SIMULATION_RUNS, rng=None):
    """Everything for ONE version transition: extract Delta E_i (Eq 1),
    run all six models for n_runs replicates each, and return the bias
    score R_i(m) (Eq 3) for the five non-random models."""
    rng = rng or np.random.default_rng()

    delta_E = observed_delta_E(G_i, G_i1)
    k = len(delta_E)
    print(f"Observed update: {k} new edges among proteins present in both versions")

    ig = IndexedGraph.build(G_i)

    jaccards = {m: [] for m in MODELS}
    for model in MODELS:
        # HCE/LCE/RWE/CRWE: precompute once per transition, reuse across
        # runs. RNBE: redraw the node attribute at the start of EVERY run.
        # RE: no per-node scores needed at all. This matches the
        # "Simulation Procedure" paragraph in Methods.
        fixed_scores = precompute_node_scores(ig, model, rng=rng) if model not in ("RE", "RNBE") else None

        for _r in range(n_runs):
            run_scores = fixed_scores
            if model == "RNBE":
                run_scores = None  # sample_new_edges draws a fresh Uniform(0,1) attribute
            sim = sample_new_edges(ig, k, model, node_scores=run_scores, rng=rng)
            jaccards[model].append(jaccard(sim, delta_E))

    R = {m: bias_score(jaccards[m], jaccards["RE"]) for m in MODELS if m != "RE"}
    return R, jaccards


def main():
    rng = np.random.default_rng(42)
    G_i, G_i1 = load_toy_transition(seed=0)
    R, jaccards = run_transition(G_i, G_i1, n_runs=N_SIMULATION_RUNS, rng=rng)

    print("\nMean Jaccard similarity per model (Eq 2):")
    for model in MODELS:
        vals = jaccards[model]
        print(f"  {model:5s}  mean J = {np.mean(vals):.4f}  (n={len(vals)} runs)")

    print("\nBias score R_i(m) relative to RE (Eq 3):")
    for model, r in sorted(R.items(), key=lambda kv: -kv[1]):
        flag = "  <- above random baseline" if r > 1 else ""
        print(f"  {model:5s}  R = {r:.3f}{flag}")


if __name__ == "__main__":
    main()
