"""
bias_score.py

Implements Eq (1)-(3) of the Methods: the empirical update Delta E_i, the
Jaccard overlap between a simulated and the real update, and the
normalized bias score R_i(m).

    Delta E_i = E_{i+1} \\ E_i                                          (Eq 1)

    J_i^(r)(m) = |Delta_tilde_E_i^(r)(m) intersect Delta E_i|
                 / |Delta_tilde_E_i^(r)(m) union Delta E_i|              (Eq 2)

    R_i(m) = <J_i^(r)(m)>_r / <J_i^(r)(RE)>_r                            (Eq 3)

IMPORTANT (matches the discussion earlier in this project about Eq 1):
Delta E_i is restricted to protein pairs that are present in BOTH G_i and
G_{i+1} (i.e. edges added between proteins that already existed in the
earlier version). New edges to a node that only appears in G_{i+1} are not
part of Delta E_i, because none of the six edge-addition models can ever
propose them -- they all sample from the complement of G_i on V_i, a fixed
node set. `observed_delta_E` enforces this restriction explicitly rather
than leaving it as an implicit side effect of set arithmetic.
"""

from __future__ import annotations

from typing import Iterable

import networkx as nx


def observed_delta_E(G_i: nx.Graph, G_i1: nx.Graph) -> set:
    """The empirical update Delta E_i (Eq 1), restricted to node pairs
    present in both G_i and G_{i+1}.

    Parameters
    ----------
    G_i, G_i1 : nx.Graph
        The earlier- and later-version networks (G_i and G_{i+1}).

    Returns
    -------
    set of frozenset({u, v}) -- edges present in G_{i+1} but not G_i,
    with both endpoints already present in V_i.
    """
    V_i = set(G_i.nodes())
    E_i = {frozenset(e) for e in G_i.edges()}
    E_i1 = {frozenset(e) for e in G_i1.edges() if set(e) <= V_i}
    return E_i1 - E_i


def jaccard(A: Iterable, B: Iterable) -> float:
    """Jaccard similarity |A ∩ B| / |A ∪ B| (Eq 2). Returns 0.0 if both
    sets are empty (rather than raising a ZeroDivisionError), since an
    empty simulated update against an empty real update carries no signal
    either way."""
    A, B = set(A), set(B)
    union = A | B
    if not union:
        return 0.0
    return len(A & B) / len(union)


def bias_score(J_model_runs: Iterable[float], J_RE_runs: Iterable[float]) -> float:
    """Normalized bias score R_i(m) = <J_i^(r)(m)>_r / <J_i^(r)(RE)>_r (Eq 3).

    Parameters
    ----------
    J_model_runs : the per-run Jaccard similarities J_i^(r)(m) for model m
        (one float per simulation replicate r).
    J_RE_runs : the per-run Jaccard similarities J_i^(r)(RE) for the same
        transition, from the SAME set of replicates used for the RE model.

    Returns
    -------
    float. R_i(m) = 1 corresponds to random expectation; > 1 indicates
    stronger agreement with the empirical update than random edge
    addition; < 1 indicates weaker agreement.
    """
    J_model_runs = list(J_model_runs)
    J_RE_runs = list(J_RE_runs)
    mean_RE = sum(J_RE_runs) / len(J_RE_runs)
    if mean_RE == 0:
        raise ValueError(
            "Mean Jaccard similarity for RE is 0 -- R_i(m) is undefined "
            "for this transition (check that k = |Delta E_i| is large "
            "enough, relative to the complement edge set, for RE to ever "
            "recover any true positives)."
        )
    mean_model = sum(J_model_runs) / len(J_model_runs)
    return mean_model / mean_RE
