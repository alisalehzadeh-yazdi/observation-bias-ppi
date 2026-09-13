"""
edge_addition_models.py

Reference implementation of the six edge-addition (interactome growth)
models described in the Methods of:

    Salehzadeh-Yazdi et al., "A quantitative analysis of observation
    biases in protein interaction networks."

Given the earlier-version network G_i = (V_i, E_i), every model assigns a
sampling weight w_m(u, v) to each candidate edge (u, v) in the complement
edge set

    Ebar_i = { (u, v) : u, v in V_i, u != v, (u, v) not in E_i }        (see Eq 4)

and draws exactly k = |Delta E_i| edges from Ebar_i, without replacement,
with probability

    p_m(u, v) = w_m(u, v) / sum_{(u', v') in Ebar_i} w_m(u', v')          (Eq 4)

Six models are implemented, matching Table 1 / Eq 5-11 of the manuscript:

    RE    Random edge addition               Eq (5)
    HCE   Highly connected edge addition     Eq (6)
    LCE   Lowly connected edge addition      Eq (7)
    RNBE  Randomized node-based edge addition Eq (8)
    RWE   Random walk edge addition          Eq (9)
    CRWE  Compact random walk edge addition  Eq (10)-(11)

Per the "Simulation Procedure" paragraph: for RE, HCE, LCE, RWE, and CRWE,
edge weights are precomputed ONCE per version transition (from G_i) and
reused across all simulation runs. For RNBE, the node attributes r_v are
redrawn independently at the start of EVERY run. `precompute_node_scores`
implements the "once per transition" step; pass `rng=None` (a fresh draw)
into `sample_new_edges` for RNBE on every run.

-----------------------------------------------------------------------
A note on exact vs. scalable sampling
-----------------------------------------------------------------------
Enumerating Ebar_i explicitly costs O(|V_i|^2), which is fine for a few
thousand nodes but not for the ~27,000-node BioGRID human network used in
this study. `sample_new_edges` therefore uses an exact-but-scalable
"generate-and-reject" scheme rather than materializing Ebar_i:

  * RE (uniform weights): draw a uniformly random pair of nodes, keep it if
    it is not a self-pair, not already an edge in G_i, and not already
    drawn in this run; repeat until k edges are collected. This is exactly
    Eq (5): every valid pair is equally likely to be accepted.

  * HCE, LCE, RWE, CRWE (weights of the additive form w(u,v) = a_u + a_v):
    a pair with probability proportional to (a_u + a_v) can be generated
    exactly, without enumerating all pairs, as follows. Draw an ordered
    pair (X, Y) by: with probability 1/2, X is drawn from a categorical
    distribution over V_i weighted by a, and Y is drawn uniformly at
    random from V_i \\ {X}; with probability 1/2, do the same with the
    roles of X and Y swapped. For P(X, Y) built this way,
    P({u, v}) = P(X=u,Y=v) + P(X=v,Y=u) = (a_u + a_v) / (A * (n - 1)),
    where A = sum_v a_v -- i.e. exactly proportional to (a_u + a_v), with
    no O(n^2) enumeration needed. The usual rejection (self-pair, existing
    edge, already drawn this run) is then applied.

  * RNBE (multiplicative weights w(u,v) = r_u * r_v): X and Y are drawn
    independently, each from the categorical distribution over V_i
    weighted by r (redrawing Y if Y == X). This is the standard
    "Chung-Lu style" construction for product-form edge weights; it is
    asymptotically exact as |V_i| grows (the small O(1/|V_i|) self-pair
    correction is handled by the resample-on-collision step, not ignored).

For small graphs (|V_i| <= EXACT_ENUMERATION_LIMIT), `sample_new_edges`
instead enumerates Ebar_i directly and samples via the Efraimidis-Spirakis
weighted-reservoir algorithm, which is the textbook-exact method for
weighted sampling without replacement. `tests/test_models.py` checks that
the two code paths agree in distribution on small graphs, so the scalable
path used for the full-size networks in the paper is validated against
the brute-force definition rather than trusted on faith.

This is a clean-room implementation of the equations as stated in the
manuscript. If you (the authors) still have the original analysis scripts
that produced Figs 2-5, treat those as the bit-exact source of truth for
the published numbers; this module is meant to be the version that lives
in the public reproducibility repository, and a clear, citable reference
that a reviewer can check line-by-line against Eq (5)-(11).
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Optional

import networkx as nx
import numpy as np

MODELS = ("RE", "HCE", "LCE", "RNBE", "RWE", "CRWE")

# Below this many nodes, sample_new_edges enumerates the complement edge
# set explicitly (used for testing / small examples). Above it, the
# scalable generate-and-reject scheme is used.
EXACT_ENUMERATION_LIMIT = 400


# ---------------------------------------------------------------------
# Node indexing helpers
# ---------------------------------------------------------------------

@dataclass
class IndexedGraph:
    """A networkx Graph plus a fixed integer index 0..n-1 for its nodes."""

    G: nx.Graph
    nodes: list
    index: dict  # node label -> integer index
    adjacency_idx: list  # adjacency_idx[i] = list of integer-index neighbors of node i
    edge_index_set: set  # set of frozenset({i, j}) for existing edges, by index

    @classmethod
    def build(cls, G: nx.Graph) -> "IndexedGraph":
        nodes = list(G.nodes())
        index = {node: i for i, node in enumerate(nodes)}
        adjacency_idx = [
            [index[nbr] for nbr in G.neighbors(node)] for node in nodes
        ]
        edge_index_set = {
            frozenset((index[u], index[v])) for u, v in G.edges()
        }
        return cls(G, nodes, index, adjacency_idx, edge_index_set)

    @property
    def n(self) -> int:
        return len(self.nodes)

    def label(self, i: int):
        return self.nodes[i]


# ---------------------------------------------------------------------
# Per-node scores for each model (Eq 6-10)
# ---------------------------------------------------------------------

def compact_random_walk_scores(
    ig: IndexedGraph,
    n_walks: int = 5000,
    walk_length: int = 100,
    p_r: float = 0.5,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Compact random walk score s_v for every node (Eq 10).

    A compact random walk of length `walk_length` starts at a node s drawn
    uniformly at random from V_i. At every step t > 0 the walker returns to
    s with probability p_r, or moves to a uniformly chosen neighbor of the
    CURRENT node with probability 1 - p_r; the walk terminates early if the
    current node has no neighbors. This is repeated for `n_walks`
    independent walks; s_v is the fraction of all recorded visits (across
    all walks, including each walk's starting position at t=0) spent at v.

    Returns
    -------
    np.ndarray of shape (n,), s_v for each node index v. Sums to 1.
    """
    rng = rng or np.random.default_rng()
    n = ig.n
    visits = np.zeros(n, dtype=np.float64)
    total_steps = 0

    for _ in range(n_walks):
        s = rng.integers(n)
        current = s
        visits[current] += 1
        total_steps += 1
        for _t in range(walk_length):
            neighbors = ig.adjacency_idx[current]
            if not neighbors:
                break  # walk terminates: current node has no neighbors
            if rng.random() < p_r:
                current = s
            else:
                current = neighbors[rng.integers(len(neighbors))]
            visits[current] += 1
            total_steps += 1

    return visits / total_steps


def precompute_node_scores(
    ig: IndexedGraph,
    model: str,
    rng: Optional[np.random.Generator] = None,
    crwe_kwargs: Optional[dict] = None,
) -> Optional[np.ndarray]:
    """Per-node score array a_v used to build w_m(u,v), for one model.

    RE needs no per-node scores (uniform weights) -> returns None.
    RNBE's attribute is redrawn per simulation run, so it is NOT computed
    here -- see `sample_new_edges`, which draws it fresh every call unless
    you pass `node_scores` explicitly (e.g. to reuse a fixed draw).
    """
    if model == "RE":
        return None
    if model == "HCE":
        return np.array([len(nbrs) for nbrs in ig.adjacency_idx], dtype=np.float64)
    if model == "LCE":
        deg = np.array([len(nbrs) for nbrs in ig.adjacency_idx], dtype=np.float64)
        k_max = deg.max() if deg.size else 0.0
        return k_max / (deg + 1.0)
    if model == "RWE":
        pr = nx.pagerank(ig.G, alpha=0.85)
        return np.array([pr[ig.label(i)] for i in range(ig.n)], dtype=np.float64)
    if model == "CRWE":
        kwargs = crwe_kwargs or {}
        return compact_random_walk_scores(ig, rng=rng, **kwargs)
    if model == "RNBE":
        rng = rng or np.random.default_rng()
        return rng.uniform(0.0, 1.0, size=ig.n)
    raise ValueError(f"Unknown model: {model!r}. Expected one of {MODELS}.")


# ---------------------------------------------------------------------
# Sampling k new edges under a given model
# ---------------------------------------------------------------------

_ADDITIVE_MODELS = {"HCE", "LCE", "RWE", "CRWE"}
_MULTIPLICATIVE_MODELS = {"RNBE"}


def _weighted_categorical_indices(weights: np.ndarray, m: int, rng: np.random.Generator) -> np.ndarray:
    """Draw m indices i.i.d. from the categorical distribution ~ weights."""
    total = weights.sum()
    if total <= 0:
        # Degenerate (e.g. all-zero weights, such as an edgeless graph for
        # HCE): fall back to uniform, matching the p_r=0 spirit of Eq 4-9
        # (a zero-weight candidate set has no informative structure left).
        return rng.integers(0, weights.size, size=m)
    cdf = np.cumsum(weights) / total
    u = rng.random(m)
    return np.searchsorted(cdf, u, side="right").clip(max=weights.size - 1)


def _generate_candidate_pairs_additive(
    weights: np.ndarray, batch: int, rng: np.random.Generator
) -> np.ndarray:
    """Batch-generate ~ (a_u + a_v) candidate pairs (see module docstring)."""
    n = weights.size
    x = _weighted_categorical_indices(weights, batch, rng)
    y = rng.integers(0, n, size=batch)
    clash = y == x
    while clash.any():
        y[clash] = rng.integers(0, n, size=clash.sum())
        clash = y == x
    swap = rng.random(batch) < 0.5
    u = np.where(swap, y, x)
    v = np.where(swap, x, y)
    return np.stack([u, v], axis=1)


def _generate_candidate_pairs_multiplicative(
    weights: np.ndarray, batch: int, rng: np.random.Generator
) -> np.ndarray:
    """Batch-generate ~ (r_u * r_v) candidate pairs (Chung-Lu style)."""
    n = weights.size
    u = _weighted_categorical_indices(weights, batch, rng)
    v = _weighted_categorical_indices(weights, batch, rng)
    clash = u == v
    while clash.any():
        v[clash] = _weighted_categorical_indices(weights, clash.sum(), rng)
        clash = u == v
    return np.stack([u, v], axis=1)


def _generate_candidate_pairs_uniform(n: int, batch: int, rng: np.random.Generator) -> np.ndarray:
    u = rng.integers(0, n, size=batch)
    v = rng.integers(0, n, size=batch)
    clash = u == v
    while clash.any():
        v[clash] = rng.integers(0, n, size=clash.sum())
        clash = u == v
    return np.stack([u, v], axis=1)


def _sample_new_edges_scalable(
    ig: IndexedGraph, k: int, model: str, node_scores: Optional[np.ndarray], rng: np.random.Generator
) -> set:
    existing = ig.edge_index_set
    chosen: set = set()
    batch = max(256, 4 * k)
    guard = 0
    while len(chosen) < k:
        guard += 1
        if guard > 10_000:
            raise RuntimeError(
                "sample_new_edges: could not find enough valid candidate "
                "edges after many attempts -- the complement edge set may "
                "be (nearly) exhausted for this k."
            )
        if model == "RE":
            pairs = _generate_candidate_pairs_uniform(ig.n, batch, rng)
        elif model in _ADDITIVE_MODELS:
            pairs = _generate_candidate_pairs_additive(node_scores, batch, rng)
        elif model in _MULTIPLICATIVE_MODELS:
            pairs = _generate_candidate_pairs_multiplicative(node_scores, batch, rng)
        else:
            raise ValueError(f"Unknown model: {model!r}. Expected one of {MODELS}.")

        for a, b in pairs:
            fs = frozenset((int(a), int(b)))
            if fs in existing or fs in chosen:
                continue
            chosen.add(fs)
            if len(chosen) == k:
                break
    return chosen


def _sample_new_edges_exact(
    ig: IndexedGraph, k: int, model: str, node_scores: Optional[np.ndarray], rng: np.random.Generator
) -> set:
    """Brute-force reference path: enumerate Ebar_i, weight it, sample
    without replacement via Efraimidis-Spirakis keys. Only used for small
    graphs (see EXACT_ENUMERATION_LIMIT) -- this is the O(n^2) method the
    scalable path above is validated against in tests/test_models.py.
    """
    n = ig.n
    candidates = [
        (i, j)
        for i, j in itertools.combinations(range(n), 2)
        if frozenset((i, j)) not in ig.edge_index_set
    ]
    if k > len(candidates):
        raise ValueError(
            f"Requested k={k} new edges but only {len(candidates)} "
            f"candidate non-edges exist."
        )
    if model == "RE":
        w = np.ones(len(candidates))
    elif model in _ADDITIVE_MODELS:
        w = np.array([node_scores[i] + node_scores[j] for i, j in candidates])
    elif model in _MULTIPLICATIVE_MODELS:
        w = np.array([node_scores[i] * node_scores[j] for i, j in candidates])
    else:
        raise ValueError(f"Unknown model: {model!r}. Expected one of {MODELS}.")

    w = np.clip(w, a_min=1e-300, a_max=None)  # keys undefined at w=0
    keys = rng.random(len(candidates)) ** (1.0 / w)  # Efraimidis-Spirakis
    top_k = np.argpartition(-keys, k - 1)[:k]
    return {frozenset(candidates[t]) for t in top_k}


def sample_new_edges(
    ig: IndexedGraph,
    k: int,
    model: str,
    node_scores: Optional[np.ndarray] = None,
    rng: Optional[np.random.Generator] = None,
) -> set:
    """Draw k new edges under `model`, from the complement of G_i (Eq 4).

    Parameters
    ----------
    ig : IndexedGraph
        The earlier-version network G_i, already indexed via
        `IndexedGraph.build`.
    k : int
        Number of edges to add (= |Delta E_i| for the real transition
        being simulated).
    model : one of MODELS
    node_scores : np.ndarray or None
        Precomputed per-node scores from `precompute_node_scores`. Ignored
        for RE. Required for HCE/LCE/RWE/CRWE (precompute once per
        transition and reuse across runs). For RNBE, pass None to draw a
        fresh Uniform(0,1) attribute for THIS run (matching the paper's
        "node attributes were drawn independently at the start of each
        run"), or pass a fixed array to reuse across calls.
    rng : numpy.random.Generator or None

    Returns
    -------
    set of frozenset({label_u, label_v}) -- the simulated new edges, using
    the ORIGINAL node labels of ig.G (not integer indices).
    """
    if model not in MODELS:
        raise ValueError(f"Unknown model: {model!r}. Expected one of {MODELS}.")
    rng = rng or np.random.default_rng()

    if model == "RNBE" and node_scores is None:
        node_scores = rng.uniform(0.0, 1.0, size=ig.n)

    if ig.n <= EXACT_ENUMERATION_LIMIT:
        idx_edges = _sample_new_edges_exact(ig, k, model, node_scores, rng)
    else:
        idx_edges = _sample_new_edges_scalable(ig, k, model, node_scores, rng)

    return {frozenset((ig.label(i), ig.label(j))) for i, j in (tuple(fs) for fs in idx_edges)}


def simulate_update(
    G_i: nx.Graph,
    k: int,
    model: str,
    ig: Optional[IndexedGraph] = None,
    node_scores: Optional[np.ndarray] = None,
    rng: Optional[np.random.Generator] = None,
) -> set:
    """Convenience one-shot wrapper: build the index (if not supplied),
    compute node scores (if not supplied), and draw k new edges.

    For repeated simulation runs on the SAME transition (as in the paper's
    n = 100 replicates), build `ig` once and precompute `node_scores` once
    for RE/HCE/LCE/RWE/CRWE, then call `sample_new_edges` directly per run
    -- recomputing PageRank or the compact random walk on every replicate
    is unnecessary work and does not match the "precomputed once per
    transition" procedure described in Methods.
    """
    ig = ig or IndexedGraph.build(G_i)
    if node_scores is None and model != "RNBE":
        node_scores = precompute_node_scores(ig, model, rng=rng)
    return sample_new_edges(ig, k, model, node_scores=node_scores, rng=rng)
