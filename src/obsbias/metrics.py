"""
Jaccard overlap and the normalized bias score, exactly as defined in the
manuscript's "Observation bias evaluation" subsection:

    J_i^(r)(m) = |dE_i_tilde^(r)(m) intersect dE_i| / |dE_i_tilde^(r)(m) union dE_i|
    R_i(m)     = <J_i^(r)(m)>_r / <J_i^(r)(RE)>_r

where the average <.>_r is taken over the `n` simulation replicates.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np


def jaccard_similarity(simulated_edges: set, empirical_edges: set) -> float:
    """J = |intersection| / |union| of two edge sets."""
    union = simulated_edges | empirical_edges
    if not union:
        return float("nan")
    return len(simulated_edges & empirical_edges) / len(union)


def bias_score(model_jaccards: Iterable[float], re_jaccards: Iterable[float]) -> float:
    """
    R_i(m) = mean(J_i(m)) / mean(J_i(RE))

    This is the quantity plotted in Figs. 3-5 of the manuscript. It is
    NOT the same as `1 / mean(J(m))`, which is what the original driver
    script computed under the name "Ratios (Actual/Model)" -- see
    docs/code_consistency_review.md section 2.1.
    """
    model_mean = float(np.mean(list(model_jaccards)))
    re_mean = float(np.mean(list(re_jaccards)))
    if re_mean == 0:
        return float("nan")
    return model_mean / re_mean
