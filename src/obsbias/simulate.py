"""
Orchestrates one database-version transition: builds the candidate edge
pool, precomputes model weights, runs `n_runs` weighted-sampling
replicates per model, and returns both the raw Jaccard values (for
plotting, cf. Fig. 2/S1) and the summary bias scores R_i(m) (cf. Figs. 3-5).
"""

from __future__ import annotations

from itertools import combinations

import networkx as nx
import numpy as np
import pandas as pd

from .growth_models import GrowthModel, sample_node_based_sequential
from .metrics import bias_score, jaccard_similarity


def _canonical_edge_set(graph: nx.Graph) -> set[tuple]:
    """
    Canonicalize every edge to `tuple(sorted((u, v)))`.

    NetworkX does not guarantee `G.edges()` returns pairs in the same
    (u, v) ordering that `itertools.combinations(nodes, 2)` would use, so
    without this step `all_possible_edges - edges` can fail to remove
    edges that already exist (see docs/code_consistency_review.md #3.1).
    """
    return {tuple(sorted(e)) for e in graph.edges()}


def build_candidate_edges(graph_early: nx.Graph) -> list[tuple]:
    """
    All non-existing pairs among graph_early's nodes: the complement edge
    set E-bar_i, returned as canonical `tuple(sorted((u, v)))` pairs.

    `itertools.combinations(nodes, 2)` yields pairs in node-list (insertion)
    order, e.g. (B, A) if "B" was added to the graph before "A" -- it does
    NOT sort each pair. Every candidate pair is therefore re-canonicalized
    here before the existing-edge check, otherwise a pair could dodge
    exclusion purely because of insertion order (the same bug class as
    #3.1 in docs/code_consistency_review.md, one level deeper).
    """
    existing = _canonical_edge_set(graph_early)
    nodes = list(graph_early.nodes())
    return [e for pair in combinations(nodes, 2) if (e := tuple(sorted(pair))) not in existing]


def get_transition_edges(graph_early: nx.Graph, graph_late: nx.Graph) -> tuple[set, set]:
    """dE_i = E_{i+1} \\ E_i (added) and E_i \\ E_{i+1} (removed), both canonicalized."""
    edges_early = _canonical_edge_set(graph_early)
    edges_late = _canonical_edge_set(graph_late)
    return edges_late - edges_early, edges_early - edges_late


def _proteins_from_edges(edges: set[tuple]) -> str:
    proteins = set()
    for u, v in edges:
        proteins.add(str(u))
        proteins.add(str(v))
    return ";".join(sorted(proteins))


def run_transition(
    graph_early: nx.Graph,
    graph_late: nx.Graph,
    models: dict[str, GrowthModel],
    n_runs: int = 100,
    include_node_based_sequential: bool = False,
    random_seed: int | None = None,
) -> dict:
    """
    Run every model in `models` for `n_runs` replicates against one
    version transition.

    Returns a dict with:
      - "jaccard": DataFrame, one column per model, one row per run
      - "proteins": DataFrame, semicolon-joined protein names in the
        intersection of predicted vs. actual added edges, per run/model
      - "bias_scores": {model_name: R_i(m)} using RE as the reference
        (requires "RE" to be one of `models`)
      - "meta": dict of transition-level bookkeeping (edge counts, etc.)

    Raises ValueError (rather than silently printing and doing nothing,
    cf. docs/code_consistency_review.md #3.3) if there are fewer candidate
    edges than the number of edges actually added.
    """
    rng = np.random.default_rng(random_seed)

    actual_added, actual_removed = get_transition_edges(graph_early, graph_late)
    candidate_edges = build_candidate_edges(graph_early)

    if len(actual_added) > len(candidate_edges):
        raise ValueError(
            f"Only {len(candidate_edges)} candidate edges available but "
            f"{len(actual_added)} edges were actually added -- cannot run "
            "a without-replacement simulation. Check that graph_early and "
            "graph_late share the intended node set."
        )

    n_edges_to_add = len(actual_added)

    # Precompute weights once per transition for every non-recomputed model
    # (matches the manuscript's Simulation Procedure for RE/HCE/LCE/RWE/CRWE).
    static_probs: dict[str, np.ndarray] = {}
    for name, model in models.items():
        if not model.recompute_each_run:
            w = model.weight_fn(graph_early, candidate_edges, rng)
            static_probs[name] = w / w.sum()

    jaccard_rows: dict[str, list[float]] = {name: [] for name in models}
    protein_rows: dict[str, list[str]] = {name: [] for name in models}
    if include_node_based_sequential:
        jaccard_rows["NBRE_sequential_experimental"] = []
        protein_rows["NBRE_sequential_experimental"] = []

    candidate_array = np.empty(len(candidate_edges), dtype=object)
    candidate_array[:] = candidate_edges

    for _ in range(n_runs):
        for name, model in models.items():
            if model.recompute_each_run:
                w = model.weight_fn(graph_early, candidate_edges, rng)
                probs = w / w.sum()
            else:
                probs = static_probs[name]

            idx = rng.choice(len(candidate_edges), size=n_edges_to_add, replace=False, p=probs)
            sampled = set(candidate_array[idx])
            intersection = sampled & actual_added

            jaccard_rows[name].append(jaccard_similarity(sampled, actual_added))
            protein_rows[name].append(_proteins_from_edges(intersection))

        if include_node_based_sequential:
            sampled = sample_node_based_sequential(graph_early, candidate_edges, n_edges_to_add, rng)
            intersection = sampled & actual_added
            jaccard_rows["NBRE_sequential_experimental"].append(jaccard_similarity(sampled, actual_added))
            protein_rows["NBRE_sequential_experimental"].append(_proteins_from_edges(intersection))

    jaccard_df = pd.DataFrame(jaccard_rows)
    jaccard_df.insert(0, "run", range(1, n_runs + 1))

    protein_df = pd.DataFrame(protein_rows)
    protein_df.insert(0, "run", range(1, n_runs + 1))

    bias_scores = {}
    if "RE" in jaccard_rows:
        for name in jaccard_rows:
            bias_scores[name] = bias_score(jaccard_rows[name], jaccard_rows["RE"])

    meta = {
        "n_nodes": graph_early.number_of_nodes(),
        "n_edges_early": graph_early.number_of_edges(),
        "n_edges_late": graph_late.number_of_edges(),
        "n_edges_added": len(actual_added),
        "n_edges_removed": len(actual_removed),
        "n_candidate_edges": len(candidate_edges),
    }

    return {"jaccard": jaccard_df, "proteins": protein_df, "bias_scores": bias_scores, "meta": meta}
