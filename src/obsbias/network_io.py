"""
Network loading and preprocessing.

Implements the preprocessing described in the manuscript's "Data sources"
paragraph:
    - self-loops and duplicate edges are removed
    - only the largest connected component (LCC) of each network version
      is retained
    - for STRING, interactions can additionally be filtered by
      `combined_score` at one of the thresholds used in the paper
      (all interactions, >400, >900, <400)

IMPORTANT: `process_file` was referenced but not defined in the original
driver script this package was cleaned from. `load_string_network` below is
a reference reimplementation reconstructed from the manuscript text for
STRING's tab/space-separated `protein1 protein2 combined_score` format only.
It has NOT been validated against the authors' original parser and should be
swapped out for it (or carefully checked against it) before being treated as
the code that produced the published figures. BioGRID, HIPPIE, and IntAct
raw files use different schemas and need their own loaders -- see the
top-level README's "Additional data/materials needed" section.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from pathlib import Path
from typing import Iterable

import networkx as nx


@dataclass(frozen=True)
class NetworkVersion:
    """One loaded, preprocessed database release."""

    path: Path
    graph: nx.Graph
    n_nodes_raw: int
    n_edges_raw: int
    n_nodes_after_lcc: int
    n_edges_after_lcc: int


def load_string_network(
    file_path: str | Path,
    score_threshold: float | None = None,
    threshold_mode: str = "greater",
    take_largest_component: bool = True,
) -> NetworkVersion:
    """
    Parse one STRING `protein.physical.links.vX.Y.txt`-style file into an
    undirected, simple graph.

    Parameters
    ----------
    file_path:
        Path to a whitespace-separated file with a header row containing at
        least `protein1`, `protein2`, and (for STRING) `combined_score`.
    score_threshold:
        If given, only keep edges passing the threshold. `None` reproduces
        the "all interactions" scheme in Fig. S1.
    threshold_mode:
        "greater" keeps `combined_score > score_threshold` (paper's `>400`,
        `>900` schemes); "less" keeps `combined_score < score_threshold`
        (paper's `<400` scheme).
    take_largest_component:
        If True (default, matches the manuscript), restrict the returned
        graph to its largest connected component.

    Returns
    -------
    NetworkVersion with both raw and post-LCC node/edge counts recorded, so
    the fraction of the network discarded by the LCC step is always
    reported (see Reviewer 1's W3 in the peer-review report: this fraction
    should be checked, not assumed negligible).
    """
    file_path = Path(file_path)
    G_raw = nx.Graph()

    with open(file_path) as fh:
        header = fh.readline().strip().split()
        try:
            p1_idx = header.index("protein1")
            p2_idx = header.index("protein2")
        except ValueError as exc:
            raise ValueError(
                f"{file_path}: expected a header with 'protein1' and "
                f"'protein2' columns, got {header!r}"
            ) from exc
        score_idx = header.index("combined_score") if "combined_score" in header else None

        for line in fh:
            parts = line.split()
            if not parts:
                continue
            u, v = parts[p1_idx], parts[p2_idx]
            if u == v:
                continue  # drop self-loops
            if score_threshold is not None and score_idx is not None:
                score = float(parts[score_idx])
                if threshold_mode == "greater" and not (score > score_threshold):
                    continue
                if threshold_mode == "less" and not (score < score_threshold):
                    continue
            G_raw.add_edge(u, v)  # networkx dedupes parallel edges automatically

    n_nodes_raw, n_edges_raw = G_raw.number_of_nodes(), G_raw.number_of_edges()

    if take_largest_component and G_raw.number_of_nodes() > 0:
        largest_cc = max(nx.connected_components(G_raw), key=len)
        G = G_raw.subgraph(largest_cc).copy()
    else:
        G = G_raw

    return NetworkVersion(
        path=file_path,
        graph=G,
        n_nodes_raw=n_nodes_raw,
        n_edges_raw=n_edges_raw,
        n_nodes_after_lcc=G.number_of_nodes(),
        n_edges_after_lcc=G.number_of_edges(),
    )


def restrict_to_common_nodes(
    versions: dict[str, NetworkVersion],
) -> dict[str, nx.Graph]:
    """
    Reproduce the ORIGINAL driver script's behavior: intersect node sets
    across every loaded version and subgraph everything down to that
    intersection.

    NOTE: this is a STRONGER restriction than "LCC of G_i alone" and is not
    described in the manuscript's Data sources paragraph. It is kept here,
    clearly named and opt-in, because the original code applied it
    unconditionally -- see docs/code_consistency_review.md section 2.3.
    Do not enable this without confirming it matches what generated the
    published figures.
    """
    node_sets = [set(v.graph.nodes()) for v in versions.values()]
    common_nodes = reduce(lambda a, b: a & b, node_sets)
    return {key: v.graph.subgraph(common_nodes) for key, v in versions.items()}


def load_versions(
    file_paths: Iterable[str | Path],
    score_threshold: float | None = None,
    threshold_mode: str = "greater",
) -> dict[str, NetworkVersion]:
    """Load several STRING releases and report basic stats for each."""
    versions = {}
    for fp in file_paths:
        nv = load_string_network(fp, score_threshold=score_threshold, threshold_mode=threshold_mode)
        pct_kept = 100 * nv.n_nodes_after_lcc / nv.n_nodes_raw if nv.n_nodes_raw else float("nan")
        print(
            f"{fp}: raw nodes={nv.n_nodes_raw}, raw edges={nv.n_edges_raw} | "
            f"after LCC: nodes={nv.n_nodes_after_lcc} ({pct_kept:.1f}% kept), "
            f"edges={nv.n_edges_after_lcc}"
        )
        versions[str(fp)] = nv
    return versions
