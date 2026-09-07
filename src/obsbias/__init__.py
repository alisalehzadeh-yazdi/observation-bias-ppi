"""
obsbias
=======

Reference implementation of the edge-addition growth models and bias-score
framework used in:

    Salehzadeh-Yazdi A, Tuncbag N, Gursoy A, Keskin O, Hutt M-T.
    "A quantitative analysis of observation biases in protein interaction
    networks."

See docs/code_consistency_review.md for a line-by-line comparison between
this code and the manuscript's Methods section, including two open
questions the authors should confirm before this is treated as the
canonical reproduction pipeline.
"""

from .metrics import jaccard_similarity, bias_score
from .growth_models import GrowthModel, build_standard_models, build_experimental_models
from .simulate import run_transition

__all__ = [
    "jaccard_similarity",
    "bias_score",
    "GrowthModel",
    "build_standard_models",
    "build_experimental_models",
    "run_transition",
]

__version__ = "0.1.0"
