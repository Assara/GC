# Triangle seed generation — proof of concept

This standalone experiment starts from **K3** at loop order 1. It does not
change the TransientGraph2 or old hairless pipelines.

For every stored graph, enumerate triples containing at least one existing
vertex and zero, one or two new vertices. New vertices receive consecutive
labels starting at `V`. Add all missing edges of that triangle. Skip triangles
already present.

| New vertices | New edges | Increase in loop order |
|---:|---:|---:|
| 0 | 1–3 | 1–3 |
| 1 | 2–3 | 1–2 |
| 2 | 3 | 1 |

Each bucket maps a canonical full-graph edge array to one unstandardized graph
and its live-component stack. The canonical key is computed without signs using
`transient_graph2_standardizer`. The stored graph and component information are
never relabelled; on duplicate keys, the first construction is retained.

The component stack is a fixed array of `uint64_t` vertex masks, of capacity
`TRIANGLE_MAX_LOOP`: K3 starts one frame, and each push costs one loop. A new
triangle must touch the latest component through at least one existing vertex.
Adding two new vertices pushes a component containing those vertices and their
attachment vertex. Adding one new vertex extends the current component. If a
triangle touches older vertices outside the latest component, pop and merge
frames until all its existing vertices belong to the live component. Popped
vertices remain part of the merged component.

This is a construction stack, not a recomputed block decomposition. Keeping only
one history per canonical graph is experimental: different histories can allow
different next triangles. Completeness must be checked again with this rule.
There is no canonical triangle-selection rule. All edges remain covered by
triangles and no isolated vertices appear. Bivalent intermediate vertices are
allowed. Cut vertices are also allowed in intermediates: attaching a triangle
at one existing vertex is necessary to reach some seeds. Seed output requires
minimum valence three and no cut vertices.

Before standardization, compute the maximum of `components(G - v) - 1` over
vertices `v`. This is a lower bound on additional loops needed to eliminate the
cut vertices: joining `k` branches without their articulation needs at least
`k - 1` additional independent cycles. Discard a candidate if this bound exceeds
the remaining loop budget. This is a necessary condition, not a guarantee that
the graph can be completed. No fixed limit on the number of cut vertices is used.

All augmentations strictly increase loop order. Process partitions in increasing
`(loop, vertices)` order, expanding each unique graph once. Release a partition
after writing and expanding it. Starting from K3, each added loop permits at most
two new vertices, so a loop cap `L` implies a vertex cap of `2L + 1`.

The previous 0–1-new-vertex restriction missed triangle-covered graphs whose
triangles must first meet at a single vertex. The two-new-vertex attachment
removes that restriction.

From the repository root, build only:

```sh
make -f GraphGeneration/Makefile.pipeline \
  triangle-seeds \
  TRIANGLE_MAX_LOOP=6
```

Build and run:

```sh
make -f GraphGeneration/Makefile.pipeline \
  run-triangle-seeds \
  TRIANGLE_MAX_LOOP=6 \
  TRIANGLE_OUTPUT=output/triangle_seeds_L6_run1
```

The output directory must be new. It contains:

- `triangles_L<L>_V<V>.g6`: distinct triangle-covered graphs surviving the loop-budget bound.
- `seeds_L<L>_V<V>.g6`: the 2-connected subset with minimum valence three.
- `counts.tsv`: candidate counts after budget pruning, unique triangle-graph counts, `seed_graphs`
  counts and elapsed seconds per processed stage. Stage times include writing
  the current partition and constructing its children.

Only reached partitions get files. Progress is printed at stage boundaries and
approximately every five seconds between parents. For initial checks, K3 should
be the sole loop-1 graph; at loop 2 the diamond and two triangles sharing one
vertex are reachable (the latter needs remaining loop budget). K4 should appear
at loop 3. The triangle and diamond must remain absent from the minimum-valence-three
seed files. K5 is reachable at loop 6.

## Second stage: splitting at fixed loop order

`tools/split_triangle_seeds.cpp` loads the minimum-valence-three `seeds_*.g6`
files for one loop order. Each stage stores `TransientGraph2` objects. Seeds
and split children are standardized with `transient_graph2_standardizer`, which
orders vertices by increasing valence, and deduplicated together.

For each graph, split every vertex of valence greater than three, scanning the
sorted valence array backwards until reaching valence three. Require both new
vertices to have valence at least three. Pass the current
loop order as **both** split loop bounds: a split adds one vertex and its connecting
edge, preserving loop order; sharing any incident edge would increase loop order
and is therefore excluded. The preserve-valence argument is zero. No additional
split filters are applied.

Process increasing vertex counts, injecting the triangle seeds at their own
vertex count. Stop at `V = 2(L - 1)`, the minimum-valence-three bound. Bivalent
triangle graphs are not inputs to this stage: these splits cannot increase the
valence of an existing bivalent vertex to three.

From the repository root, using seeds from the corrected combined run:

```sh
make -f GraphGeneration/Makefile.pipeline \
  run-triangle-splits \
  TRIANGLE_SPLIT_LOOP=6 \
  TRIANGLE_SEED_DIR=output/triangle_comparison_L6_run4/seeds \
  TRIANGLE_SPLIT_OUTPUT=output/triangle_splits_L6_run1
```

Use `triangle-splits` instead of `run-triangle-splits` to build only. The output
directory must be new. It contains `graphs_L<L>_V<V>.g6` and `counts.tsv` with
seed records, candidates, unique graphs and stage seconds. Candidates include
incoming splits and loaded seeds, before deduplication. Stage time includes
loading seeds, writing graphs and constructing the next stage. Progress is
printed at stage boundaries and approximately every five seconds between parents.
With the live-component restriction and first-history deduplication, exact
comparison through loop 9 matched all 343,683 `geng` graphs, with no missing,
extra or duplicate classes (`output/triangle_comparison_L9_run1`). This is
empirical validation through loop 9, not a completeness proof for higher orders.
That run took 4.894 seconds for seeds and 8.222 seconds for splitting, versus
15.346 seconds for geng; compilation and set comparison were excluded.

## Combined pipeline and geng comparison

Build both stages for loop orders 3 through 6:

```sh
make -f GraphGeneration/Makefile.pipeline \
  triangle-comparison TRIANGLE_MAX_LOOP=6
```

Run both stages and the comparison:

```sh
make -f GraphGeneration/Makefile.pipeline \
  run-triangle-comparison TRIANGLE_MAX_LOOP=6 \
  TRIANGLE_COMPARISON_DIR=output/triangle_comparison_L6_run1
```

The runner generates seeds once through the cap, then runs splitting for each
loop order 3 through the cap. It runs `geng -C -d3` for the same loop/vertex
partitions, targeting 2-connected graphs of minimum valence three. `labelg`
canonicalizes both outputs for exact set comparison; missing, extra and duplicate
isomorphism classes are reported separately. No pipeline graphs are filtered out
before comparison. A mismatch returns a nonzero exit status after reporting.

The fresh run directory contains seed and split outputs, `geng/` outputs,
`comparison/` canonical sets and missing/extra graph6 files, raw `logs/`,
`counts.tsv`, `performance.tsv`, and a final `summary.txt` table. Timings separate
seed generation, splitting, their combined total, and geng. They include file
writing (and live logging for the pipeline), and exclude compilation and set
comparison. These are single-run measurements; very short timings may be noisy.
