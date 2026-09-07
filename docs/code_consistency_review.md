# Code ↔ Manuscript Consistency Review

Scope: the analysis script pasted for the STRING yeast (*S. cerevisiae*, taxon 4932)
v11.5 → v12.0 transition, threshold `combined_score < 400`, compared against the
Methods section of *A quantitative analysis of observation biases in protein
interaction networks*.

This review is written so it can be dropped into the repo (`docs/`) as a record
of what was checked before the code was cleaned, and so co-authors can confirm
or correct the two open questions flagged below (‡).

## 1. Findings that match the paper (no action needed)

| Item | Paper | Code | Verdict |
|---|---|---|---|
| RE model | uniform sampling, `p_RE(u,v)=1/|Ē_i|` | `np.random.choice(..., replace=False)` (no `p`, i.e. uniform) | ✅ Match |
| HCE weight | `w = k_u + k_v` | `node_degrees.get(u,0) + node_degrees.get(v,0)` | ✅ Match |
| LCE weight | `w = k_max/(k_u+1) + k_max/(k_v+1)` | `max_degree/(k_u+1) + max_degree/(k_v+1)` | ✅ Match |
| RNBE weight | `w = r_u · r_v`, `r_v ~ U(0,1)` redrawn **every run** | `node_random_values = np.random.random(...)` redrawn **inside** the per-run loop | ✅ Match |
| RWE weight | `w = π_u + π_v`, PageRank α=0.85 | `pagerank_scores = nx.pagerank(G, alpha=0.85)`, computed **once** per transition | ✅ Match |
| CRWE algorithm | compact random walk, return prob `p_r`, `N_walks`, length `L`, score = visits/total steps | `compact_random_walk_scores(...)` — logic matches (start at random `s`, return to `s` w.p. `p_r` for `step>0`, else move to random neighbor, terminate if no neighbors, normalize by total recorded visits) | ✅ Match |
| CRWE default parameters | `N_walks=5000, L=100, p_r=0.5` | `num_walks=5000, walk_length=100, return_prob=0.5` | ✅ Match |
| "Precompute once per transition, reuse across runs" (HCE, LCE, RWE, CRWE) | stated in Simulation Procedure | `edge_probs_hcer/lcer/rwer/crwer` and `pagerank_scores`/`compact_rw_scores` are all computed **before** the `for i in range(num_runs)` loop | ✅ Match |
| Simulation replicate count | `n = 100` | `num_runs = 100` | ✅ Match |
| STRING thresholding scheme | `combined_score < 400` is one of the four schemes tested (Fig. S1) | output filename encodes `Thr_below400` | ✅ Match (this run reproduces the STRING<400 panel of Fig. S1) |

## 2. Discrepancies that need fixing before this becomes the "official" repro code

### 2.1 The bias score reported in the paper is **not** the quantity the script computes (high priority)
The paper defines everything in terms of
`R_i(m) = ⟨J_i(m)⟩_r / ⟨J_i(RE)⟩_r` (Methods, "Observation bias evaluation").
The script instead computes, per model:
```python
jaccard_actual_added = 1.0
ratio = jaccard_actual_added / mean   # i.e. 1 / mean(J(model))
```
This is not the paper's normalized bias score — it doesn't even reference the
RE model's Jaccard value. As written, running this cell would **not**
reproduce the numbers behind Figs. 2–5. This has been fixed in the cleaned
code (`metrics.bias_score`), which computes `mean(J(model)) / mean(J(RE))`.

### 2.2 Two models are simulated that are not in the manuscript's Table 1 (needs a decision, not just a fix) ‡
The paper's Table 1 lists exactly six models: RE, HCE, LCE, RNBE, RWE, CRWE.
The script computes **eight**:

| Script variable / abbreviation | What it actually does | Paper's six models? |
|---|---|---|
| `jaccard_random` / `URA` | uniform sampling | = **RE** |
| `jaccard_hcer` / `HRA` | degree-sum weighting | = **HCE** |
| `jaccard_lcer` / `PPA` | inverse-degree weighting | = **LCE** |
| `jaccard_rnber` / `LADA` | `r_u·r_v`, fresh `U(0,1)` per run | = **RNBE** |
| `jaccard_rwer` / `PWA` | PageRank-sum weighting | = **RWE** |
| `jaccard_crwer` / `CRWA` | compact-random-walk-sum weighting | = **CRWE** |
| `jaccard_node_based_random` / `SNC` | **sequential** algorithm: pick a random node, then a uniformly random *available* neighbor for that node, repeat | **not in Table 1** |
| `jaccard_cer` / `DBA` | weight `h(k) = 1/(1+|k − mean(k)|)`, i.e. favors nodes near the *mean* degree | **not in Table 1** |

Two consequences:
1. **Naming**: the script's abbreviations (`PPA, SNC, LADA, URA, DBA, PWA, HRA, CRWA`) don't correspond to the paper's names (`RE, HCE, LCE, RNBE, RWE, CRWE`) at all. Anyone trying to match a repo CSV column to a paper figure would have no way to do so without this mapping. **Fixed** in the cleaned code — models are now named exactly as in Table 1, and the two extra models are kept but clearly separated and labeled `EXPERIMENTAL_*`.
2. **Scientific question for the authors**: were `SNC` ("node-based random", not the same construction as RNBE despite the similar name) and `DBA` (mean-degree affinity) exploratory models that were tried and dropped before the final six-model design, or were they meant to appear as extra supplementary comparisons? If the latter, the manuscript's Table 1 / Supporting Information should mention them; if the former, they should probably ship in the repo as a clearly-marked `experimental/` example rather than in the main reproduction pipeline, so a reader reproducing Figs. 2–5 isn't confused by two unexplained extra bars. **This needs a co-author decision — flagged, not resolved, in the cleaned code.**

### 2.3 Node-set restriction is stricter than what the manuscript describes ‡
The manuscript's Data sources paragraph says only: *"Self-loops and duplicate
edges were removed... Only the largest connected component of each network
was retained."* Nothing there says networks are additionally restricted to
nodes common to *every* version collected.

The pasted driver code, however, does exactly that before any transition is
analyzed:
```python
common_nodes = reduce(lambda x, y: x & y, node_sets)   # intersection over ALL loaded files
filtered_networks = {fp: G.subgraph(common_nodes) for fp, G in networks.items()}
```
i.e. if six STRING versions (v9.1…v12) are loaded together, every one of them
is first cut down to only the proteins present in *all six*, and the
transition analysis (`network_v9_1`/`network_v10` in the snippet, actually
holding the v11.5/v12.0 subgraphs) runs on that intersection.

This is a materially different, stronger filter than "LCC of `G_i`" and it
directly interacts with two Methodology-review findings from the manuscript
review (candidate-edge pool construction, and how newly-appeared proteins are
handled): if a protein is absent from even one of the six STRING releases —
which will happen for proteins added or dropped over a decade of curation —
it is excluded from **every** transition's candidate pool, not just the
transitions where it's genuinely absent.

**Action needed**: confirm which of the following actually generated the
published Figs. 2–5 numbers:
(a) each transition uses only `G_i`'s own LCC, independent of other versions
(as the manuscript text implies), or
(b) each transition uses the all-versions-common-node intersection (as this
driver script does).
The cleaned code implements **both** as an explicit `restrict_to_common_nodes`
flag (default `False`, i.e. matches the manuscript text) so this is a
one-line, documented choice rather than a silent side effect — but the
correct default should be set once this is confirmed, and the manuscript's
Methods/Data sources paragraph should say explicitly which was used.

## 3. Bugs fixed during cleaning (behavior-affecting)

1. **Edge-tuple canonicalization bug.** `set(network.edges())` and
   `set(combinations(nodes, 2))` are not guaranteed to use the same
   `(u, v)` vs. `(v, u)` ordering for a given pair — NetworkX does not
   promise edges come out in the same order as a global node ordering. Where
   the two disagree on ordering, `all_possible_edges - edges_v9_1` fails to
   remove an existing edge from the candidate pool, so a handful of
   "candidate" edges could actually already exist in the graph, and would
   never appear as genuinely novel if resampled. Fixed by canonicalizing
   every edge to `tuple(sorted((u, v)))` before any set operation.
2. **Stale/incorrect log message.** The script prints
   `"...saved to 'protein_intersections_2025_2013_Mouse_...csv'"` immediately
   after saving a **yeast** STRING<400 result to a differently-named file —
   a leftover copy-paste from a prior (mouse, BioGRID 2013→2025) run. Fixed
   by deriving the printed filename from the same variable that was written.
3. **Silent no-op on insufficient candidates.** `if len(actual_added_edges) >
   len(candidate_edges): print(...)` only prints and falls through with no
   results produced; changed to raise a `ValueError` so a mis-parameterized
   run fails loudly instead of silently producing an empty results file.
4. **O(N²) rebuild every run.** The "node-based random" (`SNC`) model
   rebuilt its full non-neighbor adjacency dictionary from scratch
   (`for node in all_nodes: for other in all_nodes: ...`) on every one of
   the 100 runs. Replaced with a one-time `O(|candidate_edges|)` build from
   the already-computed `candidate_edges` list, then a cheap per-run copy.
5. **No random seed control.** Added an optional `random_seed` parameter
   threaded through `numpy` and the stdlib `random` module so a given call
   is exactly reproducible — useful for CI smoke tests and for referees who
   want to re-run a specific figure.
6. **Hardcoded, one-off file names / dates baked into the script**
   (e.g. `..._18.02.2026.csv`, variables named `network_v9_1`/`network_v10`
   that actually hold the v11.5/v12.0 networks). Replaced with names derived
   from explicit `database`, `organism`, `version_from`, `version_to`,
   `threshold` parameters, so output filenames are always self-describing
   and variable names always match their contents.

## 4. Not fixed / needs the original `process_file`

The pasted snippet calls `process_file(file_path)` but does not define it.
Per the manuscript's Data sources paragraph, this function is where
self-loop/duplicate removal and the largest-connected-component restriction
must happen (for STRING, also the `combined_score` thresholding). The cleaned
repo ships a **reference reimplementation** for STRING physical-links files
only (`network_io.load_string_network`), reconstructed from the manuscript
text, clearly commented as such. **The authors' actual `process_file` (and
equivalents for BioGRID, HIPPIE, and IntAct, which have different raw file
formats) still need to be added** — see the "Additional data/materials
needed" section of the top-level README.
