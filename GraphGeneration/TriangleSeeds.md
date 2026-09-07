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

Each bucket uses a `linear_probe_set`. Hashing and equality inspect only the used
canonical edges. Capacity is `6 * MAX_LOOP` endpoints, since `E <= 3L`.
Two standardizer variants share expansion, component forgetting, seed selection
and the splitting pipeline. The key-based variant is the default; the custom
standardizer remains an optional experiment.

- `TRIANGLE_STANDARDIZER=custom`: the component-aware standardizer
  relabels graph and masks together. Each entry stores one canonical edge array.
- `TRIANGLE_STANDARDIZER=key` (default): the older full-graph standardizer computes a key
  independently of component history. Each entry additionally keeps the original
  working edges; the component stack stays in those original labels. Duplicate
  keys retain the first working graph and stack.

The variants use separate binaries, so switching them does not reuse the wrong
build. The earlier custom standardizer permitted isomorphic duplicates among
intermediate graphs when their stack decompositions differed, even though its
final graphs matched geng through loop 9. The older key variant provides a baseline
without this history-dependent key issue.

## Component-aware standardization

`TriangleSeedStandardizer::standardize(vertices, endpoints, components)` is a
custom unsigned standardizer. It has no dependency on `Graph`, `TransientGraph2`
or `GraphStandardizer`, and does not instantiate a different implementation for
every vertex/edge dimension.

It follows the `standardize4` individualization/refinement principle:

1. Before building partitions, compare component sizes from the two ends
   inward. At the first unequal pair, discard immediately if the reversed size
   sequence is greater; otherwise start only the forward partition with its
   direction decided. A single component always starts one forward partition.
2. Only a size sequence tied with its reversal starts two partitions, one
   measured from each end. Mark the reversed attempt and copy its marker to
   descendants; the marker never enters a weight or refinement calculation.
   Split membership groups by vertex colour before comparing initial weights.
   While directions tie, apply the same refinement, weight comparison and
   individualization to both sides in a common search.
3. Once the weights distinguish the directions, drop the losing attempts. If
   only the reversed direction remains, return an empty entry immediately. If
   only the forward direction remains, continue ordinary standardization with
   no further orientation checks. `Orientation::Both` means the search weights
   never resolved the direction.
4. Materialize surviving candidate edge arrays and select the greatest array.
   There is no final reversal decision or rejection based on materialized edges.
5. Carry the selected vertex permutation to the masks, preserving their original
   stack order. For equal edge arrays, choose the relabelling with the
   lexicographically greatest ordered mask sequence. Hashing and identity use
   only edges.

A new triangle must touch the **last component**, the single live end of the
stack, through at least one existing vertex. There is no expansion across all
smallest-size ties.

The earlier fixed-round orientation-probe experiment matched all 343,683 final graphs through loop 9
(`output/triangle_comparison_L9_orientation_run1`), with seed generation at
8.762 seconds, splitting at 8.243 seconds and geng at 15.384 seconds. However,
`labelg` found 15,786 redundant intermediate records across 31 partitions.
Reversal is not the only construction-history ambiguity. For example, the
same three triangles (zero-based vertices) can be added as `012, 123, 145`,
leaving masks `{0,1,2,3}, {1,4,5}`, or as `012, 145, 123`, where reconnection
pops the stack into a single six-vertex mask. Both constructions pass the
first/last rule. That experiment produced different custom keys for their
different partitions; the older key variant avoided those duplicates. These
results describe the earlier implementation, before the joint-direction search.

Before the inward size precheck, the joint-direction search matched all 343,683 final graphs through loop 9,
with no missing, extra or duplicate final classes
(`output/triangle_comparison_L9_joint_reversal_run1`). Seed generation took
9.401 seconds and splitting 8.273 seconds. An independent `labelg` check of the
48,935 intermediate records found 16,616 duplicates across 31 partitions; the
partition counts are in `intermediate_duplicates.tsv` in that run directory.
This change therefore passes the reversal tests and final-set comparison, while
intermediate graph deduplication remains unresolved.

The standardizer tests check complementary keep/discard decisions under stack
reversal, preservation of both directions for ties, vertex-relabel invariance,
idempotence and ordered-mask relabelling. They include 100 reproducible random
triangle-addition histories, each tested in both directions with 100 vertex
permutations per direction.

With the inward size precheck and no final reversal decision, the loop-9 run
again matched all 343,683 final classes with no missing, extra or duplicate
final classes (`output/triangle_comparison_L9_size_precheck_run1`). Seed
creation took 9.349 seconds and splitting 8.303 seconds, for 17.652 seconds
combined. The preceding joint-search timing run took 9.744 seconds for seeds;
the preceding key-based run took 3.150 seconds. These are single-run timings,
not a statistical benchmark. Tests also cover size asymmetries at the second
and third mirrored component pairs, and single-component orientation handling.

## Component updates and pruning

The component stack is a fixed array of `uint64_t` vertex masks, with at most
`TRIANGLE_MAX_LOOP` entries. Adding two new vertices pushes a component containing
those vertices and the attachment vertex. Adding one new vertex extends the
live component. Reconnecting to an earlier component pops and merges the suffix
until all existing vertices of the new triangle lie in the live component.
The standardizer relabels masks but never reorders this stack.

At the key-based deduplication entry point, before
canonical-key computation or insertion, compare component sizes from both ends
inward. Discard if the first unequal mirrored pair has the larger component on
the last-end side; keep if the first-end side is larger. Continue inward only
while sizes tie. Palindromic size sequences and single-component stacks pass.
This rule does not reorder the stack or add component history to the canonical
key. The custom variant retains its own orientation handling described above.

There is no canonical triangle-selection rule. All edges remain covered by
triangles and no isolated vertices appear. Bivalent intermediate vertices are
allowed. Cut vertices are also allowed in intermediates: attaching a triangle
at one existing vertex is necessary to reach some seeds. Seed output requires
minimum valence three; cut vertices are allowed.

## Forget completed structure

The growing collection is set **A**. Expand each unique graph by adding triangles
under the attachment rules and insert the resulting graphs into A. The earlier
proposal to collect a separate set B has been replaced; there is no B wrapper,
collection, or processing step.

Whenever a graph has no bivalent vertices, forget its component structure by
replacing the ordered stack with one mask containing every vertex. Do this
before orientation rejection and canonicalization, including at the loop cap.
Every vertex then becomes available for subsequent triangle growth. This resets
construction history; it does not assert that the graph has no actual cut vertex.

Creating a new component means attaching a triangle with two new vertices at
one existing vertex. This uses the original live-component attachment rule,
regardless of how many bivalent vertices remain. The proposed restriction to
zero bivalent vertices, or attachment at the sole bivalent vertex, was removed.
All triangle-grown vertices already have degree at least two, so "no bivalent
vertices" is equivalent to minimum degree at least three here.

Cut-budget pruning remains removed. Completed graphs with actual cut vertices
are retained in A and in the minimum-valence-three seed output. The splitting
stage resolves cut vertices before ordinary splitting.
The earlier benchmark counts and timings below predate this growth policy.

### Budget for a cubic completion

Before standardization, reject a graph with `b` bivalent vertices at loop order
`L` if `b > 2 * (MAX_LOOP - L)`. A triangle addition can remove at most two
bivalent vertices per added loop; introducing new vertices cannot improve that
bound. Vertex splits preserve loop order and do not repair bivalent vertices.
Higher degrees are allowed in intermediates because they can be split later.
A single triangle can repair three bivalent vertices, but doing so requires at
least two new edges between existing vertices and therefore at least two added
loops. The bound is per loop, not per triangle.

Also reject `V > 2 * (MAX_LOOP - 1)`: a connected cubic graph has
`V = 2 * (L - 1)`, and triangle additions and splits never decrease vertex count.
These are necessary bounds, not a guarantee of completion. They do not restore
the old cut-vertex bound or the removed attachment restriction. Existing output
still includes the stages of minimum degree three before the final cubic stage.

With this completion budget, `output/triangle_comparison_L9_completion_budget_run1`
matches all 343,683 geng graphs through loop 9, with no missing graphs, extras,
or duplicates. Seed generation takes 0.110 seconds, splitting 6.795 seconds,
and the full pipeline 6.905 seconds versus geng's 15.424 seconds (compilation
and comparison excluded). The component/budget tests pass. This replaces the
unpruned loop-9 attempt, which failed with `std::bad_alloc` after 53.631 seconds.

The loop-10 run (`output/triangle_comparison_L10_completion_budget_run1`) takes
0.616 seconds for seeds through loop 10, 215.168 seconds for loop-10 splitting
alone, and 222.481 seconds for the full pipeline through loop 10. Seed generation
previously took 32.026 seconds. Exact comparison reuses the completed canonical
geng reference from `output/triangle_comparison_L10_plain_graph_run1`; geng was
not timed again. Loop 10 produces 7,833,281 of 7,833,332 reference graphs, with
the same 51 missing classes as before and no extras or duplicates. All missing
files were checked byte-for-byte against the previous run. Through loop 10 the
totals are 8,176,964 versus 8,177,015. The new pruning introduces no additional
missing classes in that comparison, which still used the old cut-free seed selection.

With component forgetting retained and the attachment restriction removed,
`output/triangle_comparison_L6_forget_unrestricted_run1` matches all 117 geng
graphs through loop 6, with no missing graphs, extras, or duplicates. The
component-stack tests also pass.

The removed attachment restriction (`output/triangle_comparison_L6_forget_components_run1`)
matches geng at loops 3 and 4, but misses 3 final graphs at loop 5 and 6 at loop 6
(108 versus 117 total; no extras or duplicates). The smallest missing graph is
graph6 `E` followed by a backtick and `~o`: two nonadjacent vertices each joined
to four other vertices, with a matching on those four vertices. Its degrees are
`3,3,3,3,4,4` and all ten edges lie in triangles.

Its only triangles are `014`, `015`, `234`, and `235`. Starting with any of
them, a legal next addition can only complete a diamond. That diamond has two
bivalent vertices; either remaining triangle would start a new component while
leaving one old vertex bivalent, so both are forbidden. An exhaustive enumeration
of triangle additions within this target, without stack or orientation
restrictions, confirms that it is unreachable under the new attachment rule.
It also cannot come from the existing splitting stage: the new edge of a vertex
split has no common neighbor at its endpoints, whereas every edge of this target
lies in a triangle. This counterexample motivated removing the attachment restriction.

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

- `triangles_L<L>_V<V>.g6`: distinct triangle-covered graphs within the loop cap.
- `seeds_L<L>_V<V>.g6`: all graphs of minimum valence three, including cut vertices.
- `counts.tsv`: candidate counts, unique triangle-graph counts, `seed_graphs`
  counts and elapsed seconds per processed stage. Stage times include writing
  the current partition and constructing its children.

Only reached partitions get files. Progress is printed at stage boundaries and
approximately every five seconds between parents. For initial checks, K3 should
be the sole loop-1 graph; at loop 2 the diamond and two triangles sharing one
vertex are reachable. K4 should appear
at loop 3. The triangle and diamond must remain absent from the minimum-valence-three
seed files. K5 is reachable at loop 6.

## Second stage: splitting at fixed loop order

Seeds are filtered only by minimum degree three. Before expanding each graph,
`CutVertexSplitRule<V>` finds a cut vertex, if any, and records its neighbours
in each connected component of G-v. If there is a cut vertex, split only that
vertex and retain a partition only when both replacement vertices meet every
component of G-v. This removes the selected articulation without introducing
one at either replacement. Complementary partitions are still identified by
fixing the first incidence on the old vertex. The operation uses plain
`Graph::splitGraph`, may break triangles, and preserves loop order.

Recheck each child: another cut vertex must be resolved before any ordinary
vertex split. Once no cut vertices remain, use the existing splitting rules.
Only graphs without cut vertices are written to `graphs_*.g6`; pending repair
graphs remain internal. The `graphs` column in splitting counts counts written
graphs, while progress also reports `pending_cut_graphs`.

Validation: `output/triangle_comparison_L10_cut_splits_run1` matches all
7,833,332 loop-10 graphs and all 8,177,015 graphs through loop 10 against the
cached geng reference, with no missing graphs, extras, or duplicates. All 51
previously missing classes are recovered. Seeds take 0.637 seconds, splitting
through loop 10 takes 231.848 seconds, and the full pipeline takes 232.485
seconds (compilation and exact comparison excluded). Focused tests exhaust the
complement-distinct splits of two K4s sharing a vertex against direct connectivity
checks and verify that repairing one cut in a three-K4 chain leaves the next cut
for subsequent processing.

Seed-only scaling runs with the same rules (no splitting or geng comparison):

| Loop cap | Seeds at cap | Seeds through cap | Generation time | Peak RSS |
|---:|---:|---:|---:|---:|
| 11 | 8,636 | 10,821 | 4.084 s | 36.8 MiB |
| 12 | 47,750 | 58,571 | 28.029 s | 237.3 MiB |

Results are in `output/triangle_seeds_L11_completion_budget_run1` and
`output/triangle_seeds_L12_completion_budget_run1`, with sibling `.log` and
`.metrics.json` files. Times include generation through the cap and output,
but exclude compilation. Peak RSS is the generator child process's Linux
`ru_maxrss`, measured through Python's `resource` module. Both runs completed;
seed file record counts were checked against their manifests. These counts
have not been independently verified against an exhaustive reference.


`tools/split_triangle_seeds.cpp` loads the minimum-valence-three `seeds_*.g6`
files for one loop order. Each stage stores plain `Graph` objects. Seeds and
split children are standardized directly with
`GraphStandardizer::standardize_no_sign` and deduplicated in a `std::set`.
The splitting path uses neither `TransientGraph2` nor `BasisElement`/`LinComb`
operations.

For each graph, scan all vertices and split every vertex of valence greater
than three. Enumerate subsets of incident edges to move to the new vertex,
leaving incidence zero on the old vertex to identify complementary splits.
Each side receives at least two original incidences and the connecting edge,
so both child vertices have valence at least three. Call `Graph::splitGraph`
to construct each plain child graph. Every split adds one vertex and one edge,
preserving loop order. There are no shared incidences, valence-sorting pass or
additional split filters.

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
Exact comparison with the ordered stack, first/last check and one live component
matched all 343,683 graphs through loop 9, with no missing, extra or duplicate
final classes (`output/triangle_comparison_L9_ordered_stack_run1`). Seed generation
took 9.377 seconds and splitting 8.283 seconds; geng took 15.705 seconds.
Compilation and comparison are excluded. This is empirical validation through
loop 9, not a completeness proof for higher orders.

The older key variant with the same first/last rule also matched all 343,683
final graphs (`output/triangle_comparison_L9_key_first_last_run1`). Its seed stage
took 3.215 seconds, splitting 8.289 seconds, and geng 15.391 seconds. All its
intermediate triangle partitions were independently checked with `labelg` and
contained no duplicate isomorphism classes.

The default key-based variant with full inward size comparison matched all
343,683 final graphs through loop 9, with no missing, extra or duplicate final
classes (`output/triangle_comparison_L9_key_inward_run1`). Seed generation took
3.144 seconds and splitting 8.267 seconds (11.411 seconds combined). Independent
`labelg` checks found 0 duplicate classes among 28,289 intermediate records.

Direct `Graph` splitting matched all 343,683 final graphs through loop 9, with
no missing, extra or duplicate final classes
(`output/triangle_comparison_L9_direct_graph_splits_run1`). Every stage also
matched the transient baseline's candidate and distinct-graph counts. Splitting
took 6.689 seconds, compared with 8.267 seconds for the preceding unseeded
inward-check transient run (about 19% less time). Seed generation took 3.202
seconds, making the combined pipeline 9.891 seconds. These are single-run
measurements with the same `-O1` flags; compilation and comparison are excluded.

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

To run the default key-based variant with the full inward size check:

```sh
make -f GraphGeneration/Makefile.pipeline \
  run-triangle-comparison TRIANGLE_MAX_LOOP=9 \
  TRIANGLE_STANDARDIZER=key \
  TRIANGLE_COMPARISON_DIR=output/triangle_key_first_last_L9_run2
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
