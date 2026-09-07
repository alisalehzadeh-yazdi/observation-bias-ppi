"""
Edge-addition growth models.

The six models below (`RE, HCE, LCE, RNBE, RWE, CRWE`) are exactly Table 1 of
the manuscript. Two additional models found in the original driver script
(`NBRE_SEQUENTIAL`, `DCA`) are NOT in Table 1 and are kept separately under
`build_experimental_models()` -- see docs/code_consistency_review.md section
2.2. Do not include them in a run that is meant to reproduce Figs. 2-5
without confirming with the authors whether they belong in the paper's
Supporting Information.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable

import networkx as nx
import numpy as np


@dataclass
class GrowthModel:
    """
    A weighted edge-addition model.

    `weight_fn(graph, candidate_edges, rng) -> np.ndarray` returns one
    non-negative weight per entry of `candidate_edges` (same order). The
    weights are normalized to probabilities by the caller.

    `recompute_each_run` mirrors the manuscript's Simulation Procedure:
    False for RE/HCE/LCE/RWE/CRWE (edge weights "precomputed once per
    transition and reused across all simulation runs"), True for RNBE
    (node attributes "drawn independently at the start of each run").
    """

    name: str
    description: str
    weight_fn: Callable[[nx.Graph, list[tuple], np.random.Generator], np.ndarray]
    recompute_each_run: bool = False


def _weights_re(graph, candidate_edges, rng) -> np.ndarray:
    """RE: unbiased, uniform weight over the complement edge set."""
    return np.ones(len(candidate_edges), dtype=float)


def _weights_hce(graph, candidate_edges, rng) -> np.ndarray:
    """HCE: w(u, v) = k_u + k_v."""
    degree = dict(graph.degree())
    return np.array([degree.get(u, 0) + degree.get(v, 0) for u, v in candidate_edges], dtype=float)


def _weights_lce(graph, candidate_edges, rng) -> np.ndarray:
    """LCE: w(u, v) = k_max/(k_u+1) + k_max/(k_v+1)."""
    degree = dict(graph.degree())
    k_max = max(degree.values()) if degree else 0
    return np.array(
        [k_max / (degree.get(u, 0) + 1) + k_max / (degree.get(v, 0) + 1) for u, v in candidate_edges],
        dtype=float,
    )


def _weights_rnbe(graph, candidate_edges, rng) -> np.ndarray:
    """RNBE: w(u, v) = r_u * r_v, r_v ~ Uniform(0, 1), redrawn every run."""
    nodes = list(graph.nodes())
    r = dict(zip(nodes, rng.random(len(nodes))))
    return np.array([r.get(u, 0.0) * r.get(v, 0.0) for u, v in candidate_edges], dtype=float)


def _weights_rwe(graph, candidate_edges, rng) -> np.ndarray:
    """RWE: w(u, v) = pi_u + pi_v, PageRank with damping factor alpha=0.85."""
    pagerank = nx.pagerank(graph, alpha=0.85)
    return np.array([pagerank.get(u, 0.0) + pagerank.get(v, 0.0) for u, v in candidate_edges], dtype=float)


def compact_random_walk_scores(
    graph: nx.Graph,
    num_walks: int = 5000,
    walk_length: int = 100,
    return_prob: float = 0.5,
    rng: np.random.Generator | None = None,
) -> dict:
    """
    s_v = (visits to v summed over all walks) / (total recorded steps
    summed over all walks), where each walk starts at a random node s,
    at each step > 0 returns to s with probability `return_prob`, else
    moves to a uniformly random neighbor, and terminates early if the
    current node has no neighbors. Matches the manuscript's CRWE
    description and default parameters (N_walks=5000, L=100, p_r=0.5).
    """
    rng = rng or np.random.default_rng()
    nodes = list(graph.nodes())
    visit_counts: dict = defaultdict(int)

    for _ in range(num_walks):
        start = nodes[rng.integers(len(nodes))]
        current = start
        for step in range(walk_length):
            visit_counts[current] += 1
            if step > 0 and rng.random() < return_prob:
                current = start
                continue
            neighbors = list(graph.neighbors(current))
            if not neighbors:
                break
            current = neighbors[rng.integers(len(neighbors))]

    total_visits = sum(visit_counts.values())
    scores = {node: visit_counts.get(node, 0) / total_visits for node in nodes}
    return scores


def _weights_crwe_factory(num_walks: int, walk_length: int, return_prob: float):
    def _weights_crwe(graph, candidate_edges, rng) -> np.ndarray:
        scores = compact_random_walk_scores(
            graph, num_walks=num_walks, walk_length=walk_length, return_prob=return_prob, rng=rng
        )
        return np.array([scores.get(u, 0.0) + scores.get(v, 0.0) for u, v in candidate_edges], dtype=float)

    return _weights_crwe


def build_standard_models(
    crwe_num_walks: int = 5000, crwe_walk_length: int = 100, crwe_return_prob: float = 0.5
) -> dict[str, GrowthModel]:
    """The six models from the manuscript's Table 1, keyed by the paper's own abbreviations."""
    return {
        "RE": GrowthModel("RE", "Random edge addition (unbiased baseline)", _weights_re),
        "HCE": GrowthModel("HCE", "Highly connected edge addition", _weights_hce),
        "LCE": GrowthModel("LCE", "Lowly connected edge addition", _weights_lce),
        "RNBE": GrowthModel(
            "RNBE", "Randomized node-based edge addition", _weights_rnbe, recompute_each_run=True
        ),
        "RWE": GrowthModel("RWE", "Random walk edge addition (PageRank, alpha=0.85)", _weights_rwe),
        "CRWE": GrowthModel(
            "CRWE",
            "Compact random walk edge addition",
            _weights_crwe_factory(crwe_num_walks, crwe_walk_length, crwe_return_prob),
        ),
    }


# ---------------------------------------------------------------------------
# Experimental models present in the original driver script but NOT in the
# manuscript's Table 1. See docs/code_consistency_review.md section 2.2.
# Kept separate on purpose -- do not merge into build_standard_models()
# without a co-author decision on whether these belong in the paper.
# ---------------------------------------------------------------------------


def _weights_dca(graph, candidate_edges, rng) -> np.ndarray:
    """
    EXPERIMENTAL, not in Table 1. Originally named "CER" in the driver
    script (`DBA` in its output columns). Favors edges whose endpoints are
    close to the network's *mean* degree: h(k) = 1 / (1 + |k - mean(k)|).
    """
    degree = dict(graph.degree())
    mean_degree = np.mean(list(degree.values())) if degree else 0.0
    h = {k: 1.0 / (1.0 + abs(k - mean_degree)) for k in set(degree.values())}
    return np.array([h[degree.get(u, 0)] + h[degree.get(v, 0)] for u, v in candidate_edges], dtype=float)


def sample_node_based_sequential(
    graph: nx.Graph, candidate_edges: list[tuple], n_edges: int, rng: np.random.Generator
) -> set[tuple]:
    """
    EXPERIMENTAL, not in Table 1. Originally named "node-based random"
    model in the driver script (`SNC` in its output columns). Distinct
    from RNBE: repeatedly picks a uniformly random node that still has an
    available (non-neighbor) partner, then a uniformly random available
    partner for that node, and adds that edge -- rather than a single
    weighted draw over the whole candidate-edge set. Built from
    `candidate_edges` in O(|candidate_edges|) rather than the original
    script's O(|V|^2) nested loop.
    """
    available: dict = defaultdict(set)
    for u, v in candidate_edges:
        available[u].add(v)
        available[v].add(u)

    nodes_with_options = [n for n, opts in available.items() if opts]
    chosen: set[tuple] = set()

    while len(chosen) < n_edges and nodes_with_options:
        node = nodes_with_options[rng.integers(len(nodes_with_options))]
        options = available[node]
        if not options:
            nodes_with_options.remove(node)
            continue
        partner = list(options)[rng.integers(len(options))]
        edge = tuple(sorted((node, partner)))
        if edge in chosen:
            available[node].discard(partner)
            available[partner].discard(node)
            continue
        chosen.add(edge)
        available[node].discard(partner)
        available[partner].discard(node)
        if not available[node] and node in nodes_with_options:
            nodes_with_options.remove(node)
        if not available[partner] and partner in nodes_with_options:
            nodes_with_options.remove(partner)

    return chosen


def build_experimental_models() -> dict[str, GrowthModel]:
    """Models present in the original script but absent from the manuscript's Table 1."""
    return {
        "DCA_experimental": GrowthModel(
            "DCA_experimental",
            "Mean-degree affinity (NOT in manuscript Table 1; was 'CER'/'DBA')",
            _weights_dca,
        ),
        # NBRE_SEQUENTIAL is not weight-based, so it is not returned here;
        # call sample_node_based_sequential(...) directly (see simulate.py).
    }
