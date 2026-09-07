#!/usr/bin/env python
"""
Example / reproduction script: run the six-model observation-bias analysis
for one database version transition.

This replaces the original ad hoc notebook cell that hardcoded a single
yeast STRING v11.5 -> v12.0, `combined_score < 400` run with variable names
left over from an earlier v9.1 -> v10 test. Everything that was hardcoded
before is now a CLI argument, and output filenames are always derived from
those arguments (no more stray dates or mismatched print statements).

Example (reproduces the original snippet's STRING<400, yeast, v11.5->v12.0 run):

    python scripts/run_single_transition.py \\
        --early data/raw/string/Y4932.protein.physical.links.v11.5.txt \\
        --late  data/raw/string/Y4932.protein.physical.links.v12.0.txt \\
        --database STRING --organism yeast \\
        --version-from 11.5 --version-to 12.0 \\
        --threshold 400 --threshold-mode less \\
        --n-runs 100 --include-experimental \\
        --outdir results/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from obsbias.growth_models import build_experimental_models, build_standard_models  # noqa: E402
from obsbias.network_io import load_string_network  # noqa: E402
from obsbias.simulate import run_transition  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--early", required=True, help="Path to the earlier version's raw file")
    p.add_argument("--late", required=True, help="Path to the later version's raw file")
    p.add_argument("--database", required=True, help="e.g. STRING, BioGRID, HIPPIE, IntAct")
    p.add_argument("--organism", required=True, help="e.g. yeast, human, ecoli, drosophila, mouse")
    p.add_argument("--version-from", required=True, help="e.g. 11.5")
    p.add_argument("--version-to", required=True, help="e.g. 12.0")
    p.add_argument("--threshold", type=float, default=None, help="STRING combined_score threshold")
    p.add_argument("--threshold-mode", choices=["greater", "less"], default="greater")
    p.add_argument("--n-runs", type=int, default=100)
    p.add_argument("--random-seed", type=int, default=None)
    p.add_argument(
        "--include-experimental",
        action="store_true",
        help="Also run the two non-Table-1 models (DCA, node-based-sequential). "
        "Off by default -- see docs/code_consistency_review.md section 2.2.",
    )
    p.add_argument("--outdir", default="results")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    nv_early = load_string_network(args.early, score_threshold=args.threshold, threshold_mode=args.threshold_mode)
    nv_late = load_string_network(args.late, score_threshold=args.threshold, threshold_mode=args.threshold_mode)

    models = build_standard_models()
    if args.include_experimental:
        models.update(build_experimental_models())

    result = run_transition(
        nv_early.graph,
        nv_late.graph,
        models,
        n_runs=args.n_runs,
        include_node_based_sequential=args.include_experimental,
        random_seed=args.random_seed,
    )

    thr_tag = f"_thr{args.threshold_mode}{args.threshold}" if args.threshold is not None else ""
    tag = f"{args.database}_{args.organism}_v{args.version_from}-v{args.version_to}{thr_tag}_n{args.n_runs}"

    jaccard_path = outdir / f"jaccard_{tag}.csv"
    protein_path = outdir / f"protein_intersections_{tag}.csv"
    result["jaccard"].to_csv(jaccard_path, index=False)
    result["proteins"].to_csv(protein_path, index=False)

    print(f"Transition: {args.database} {args.organism} v{args.version_from} -> v{args.version_to}")
    print(result["meta"])
    print("\nBias scores R_i(m) (mean(J(model)) / mean(J(RE))):")
    for name, score in sorted(result["bias_scores"].items(), key=lambda kv: -kv[1]):
        print(f"  {name:>28s}: {score:.3f}")
    print(f"\nSaved: {jaccard_path}")
    print(f"Saved: {protein_path}")


if __name__ == "__main__":
    main()
