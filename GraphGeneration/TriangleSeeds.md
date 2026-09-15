# Triangle seed generation — proof of concept

This standalone experiment starts from **K3** at loop order 1. It does not
change the TransientGraph2 or old hairless pipelines.

## Separate pipeline grouped by valence array

`tools/split_triangle_seeds_by_valence.cpp` builds one target degree sequence at
a time. The existing shared-bucket splitter remains available; both use the
same graph6 I/O and bucket implementation in `TriangleSplitCommon.hpp`.

```sh
make -f GraphGeneration/Makefile.pipeline triangle-splits-valence TRIANGLE_SPLIT_LOOP=9
OMP_NUM_THREADS=8 ./build/triangle_splits_valence_L9 \
  output/triangle_comparison_L10_cut_splits_run1/seeds \
  output/triangle_splits_valence_L9_new_run
```

The command requires a fresh output directory. It writes the usual per-vertex
`graphs_L9_V*.g6` and `counts.tsv`, plus `valence_counts.tsv` with parent-group,
candidate, stored-graph, output-graph, and time counts for each target sequence.
Graph6 representatives have nondecreasing vertex degrees. Output line order can
vary with worker scheduling; graph identities are deterministic.

`ValenceGraphStandardizer.hpp` canonicalizes with the existing standardizer,
then orders that canonical graph by (degree, canonical vertex index). The
existing entry point named `standardize_no_sign_increasing_valency` does not
always produce sorted degrees, so this pipeline explicitly enforces the order.

Before reading parents, `ValenceSplitPlan.hpp` builds the parent-to-target
routing table. Splitting degree d and moving k incidences replaces d with
(k+1, d-k+1). Each parent/target pair therefore supplies the vertex indices and
moved sizes for every graph in that parent group. Both unequal moved sizes must
remain: incidence zero stays on the old vertex, so either side of the degree
pair can contain it. Only adjacency-dependent cut and reduction checks remain
per graph. These are the existing filters, including internal cut-vertex graphs.

Parent groups live in private binary files under `working/`, containing only
canonical half-edge bytes. One target hash set and a batch of at most 4096
parents occupy RAM. A completed target is drained directly to its binary file
and public graph6 output without allocating a packed copy. Earlier parent files
are deleted only after all their target groups have been processed. All scratch
files are removed on successful completion; partial runs retain their files and
cannot be resumed or overwritten by this command. Scratch disk space can include
two vertex stages in binary, in addition to the saved graph6 outputs.

The largest individual valence group must still fit in memory. Hash tables grow
normally; this is not a hard memory cap or an out-of-core hash set. No L12 run
has been attempted with this pipeline.

L9 validation and benchmark (8 workers, default -O1, three timed runs, compilation
excluded) are saved in `output/triangle_splits_L9_valence_benchmark`:

| Pipeline | Median seconds | Median peak RSS |
|---|---:|---:|
| Existing shared buckets | 0.560 | 16.68 MiB |
| One target valence group | 0.778 | 11.78 MiB |

All 326,192 output graph identities and 1,095,580 candidates match the saved
reference; one-worker and eight-worker outputs also agree. The dedicated test
checks 25,413 actual split routings, sorted degrees, canonical idempotence and
invariance under vertex relabelling. The existing pipeline's extracted common
code also reproduces its previous output. L9 shows about 29% lower RSS and 39%
longer runtime; it does not establish an order-of-magnitude memory saving.

```sh
make -f GraphGeneration/Makefile.pipeline build/valence_split_pipeline_test_L9
./build/valence_split_pipeline_test_L9 \
  output/triangle_comparison_L10_cut_splits_run1/splits_L9 \
  output/triangle_splits_L9_valence_benchmark/run2
```

### L10 valence-group benchmark

Three alternating runs with eight workers and -O1, excluding compilation and
verification, including graph-file writes:

| Pipeline | Median seconds | Median peak RSS |
|---|---:|---:|
| Existing shared buckets | 16.976 | 285.73 MiB |
| One target valence group | 26.752 | 125.48 MiB |

This is 56.1% lower peak RSS and 57.6% longer runtime. The full comparison of
valence1 against the saved L10 reference verifies all 7,833,332 graph identities,
with no duplicates, plus 50,413 split routings. All six timed runs match graph
counts at every vertex stage. Candidate totals differ by two (25,775,864 versus
25,775,862): sorted-degree labelling changes some intermediate cut-first paths,
without changing any public output graph. Results, graphs, binary hashes, and
verification logs are in `output/triangle_splits_L10_valence_benchmark`.

## Independent experiment without a component stack

`UnrestrictedTriangleSeedGenerator<MaxLoop>` is a separate class in
`GraphGeneration/UnrestrictedTriangleSeedGenerator.hpp`. The original generator
and its default build target remain unchanged. The side class stores only
canonical graph edges in `TriangleGraphEntry`, and expands from that labelling.
It has no component masks, original-edge payload, live-end restriction, stack
reversal test, or component forgetting. It retains the existing full-graph key
standardizer and cubic-completion budget. Every triangle augmentation sharing
at least one existing vertex is considered, including those with two new vertices.
Output uses the existing seed/count formats and includes cut-vertex seeds.

Build and run, for example:

```sh
make -f GraphGeneration/Makefile.pipeline unrestricted-triangle-seeds TRIANGLE_MAX_LOOP=12
./build/triangle_seeds_unrestricted_L12 output/triangle_seeds_L12_unrestricted_new_run
```

Measured seed-only runs (compilation excluded; time is through the loop cap):

| Loop cap | Stack time | No-stack time | Stack peak RSS | No-stack peak RSS | Seeds at cap, stack | Seeds at cap, no stack |
|---:|---:|---:|---:|---:|---:|---:|
| 11 | 4.084 s | 6.362 s | 36.8 MiB | 17.8 MiB | 8,636 | 8,650 |
| 12 | 28.029 s | 57.133 s | 237.3 MiB | 80.2 MiB | 47,750 | 47,819 |

Results: `output/triangle_seeds_L11_unrestricted_run1` and
`output/triangle_seeds_L12_unrestricted_run1`, with sibling logs and metrics.
Their `baseline_comparison.json` files compare canonical seed sets against the
previous component-stack runs. No baseline seeds are lost. The side class adds
19 seeds through loop 11 and 88 through loop 12: 1 at loop 9, 4 at loop 10,
14 at loop 11, and 69 at loop 12. Manifest counts match the seed files. This
demonstrates that the stack suppresses seed classes; it does not establish that
these additional seeds produce previously missing final graphs. Splitting was
not rerun for this seed-only experiment.

The remainder describes the original component-stack generator and its history.

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

### Shared buckets and memory

Stages now use 64 `linear_probe_set` buckets. OpenMP workers share them rather
than retaining duplicate children in per-worker sets. Resizing duplicates only
one bucket, and conversion to a sorted vector releases buckets as they are
copied. The serial build uses the same layout without locks.

Loop-10 benchmark, eight workers, identical seeds, split filter and `-O1` flags;
three alternating runs per variant, compilation excluded:

| Storage | Median wall time | Median peak RSS |
|---|---:|---:|
| Worker-local sets | 33.04 s | 721.0 MiB |
| Shared buckets | 22.20 s | 301.4 MiB |

All 7,833,332 graphs and all candidate counts match byte-for-byte on every run.
Raw measurements, logs and source snapshots are in
`output/triangle_splits_L10_memory_benchmark`. This measures about 58% less peak
memory and 33% less elapsed time. The loop-11 retry was interrupted when loop 10
was selected as the benchmark; loop-11 completion has not been established.

### Reserving the next stage

`SplitStageEstimate.hpp` estimates the next stage from loop order `L`, current
vertex count `V`, and the current number of unique stored graphs `N`:

`estimated_next = ceil(N * model(L, V+1) / model(L, V))`.

The model sums degree sequences of minimum degree three, weighted by inverse
degree factorials, counts endpoint pairings, and applies a sparse simplicity
correction. It ignores exact automorphism and cut-vertex effects, so this is a
heuristic sizing hint, not a bound. Degree surplus is `2*L - 2 - V`; coefficients
are computed with a small dynamic program once per stage. For stages below
1,024 graphs, use the current size as a conservative reservation instead of
extrapolating the unreliable dense/symmetric regime.

`SplitStageEstimator` now calibrates the model using the preceding stage's
observed error. Multiply the raw next estimate by
`sqrt(clamp(actual_current / previous_raw_estimate, 0.5, 2))`. This corrects half
the multiplicative error, limits the correction to approximately 0.71–1.41, and
does not feed corrected predictions back into the model. Calibration starts only
after a model-based prediction from a stage of at least 1,024 graphs is available.
State resets at the beginning of each pipeline run; no future stage counts or
loop-specific lookup tables are used.

Reserve `ceil(estimated_next / 64)` elements in each bucket before expansion.
`linear_probe_set` still rounds to power-of-two slot capacity with an 87.5%
load limit and grows normally if needed. Its insertion now preserves reserved
capacity; previously it could shrink a reservation on the next insertion.
Tests cover reservation retention, contents and growth beyond an underestimate.
This avoids repeated growth; it does not eliminate power-of-two spare capacity.

For example, the observed loop-10 V=13 size of 1,632,970 predicts 1,896,135
entries at V=14 (observed 2,046,783); V=14 predicts 1,600,047 at V=15
(observed 1,695,129). Both examples reserve the same power-of-two capacities
ultimately needed in those stages, despite underestimating graph counts.

Three alternating loop-9 runs with eight shared-bucket workers and `-O1`
gave medians of 0.579 seconds without reservation and 0.563 seconds with it.
All graph bytes and candidate counts match. Prediction/actual stage counts and
measurements are in `output/triangle_splits_L9_reserve_benchmark`. Memory at this
small scale was approximately unchanged. Loop-10 examples above evaluate the
model on saved sizes; they are not a new loop-10 timing or memory benchmark.

Replaying the saved stages through loop 10 reduces count-weighted absolute
prediction error on model-sized stages from 10.37% to 4.40% at loop 9 and from
8.39% to 2.91% at loop 10. For loop 10, predictions for V=12..16 are 834,056,
1,586,888, 1,994,314, 1,662,395 and 878,710; actual counts are 868,719,
1,632,970, 2,046,783, 1,695,129 and 900,109. These are retrospective checks,
not error guarantees for larger loops.

Three paired loop-9 shared-bucket runs retain identical graph bytes and counts.
Median times are 0.568 seconds for the raw model and 0.569 for calibration, with
approximately unchanged peak memory. The improvement is prediction accuracy;
power-of-two capacities often remain the same. Logs, predictions and replay
results are in `output/triangle_splits_L9_adaptive_estimate_benchmark`.

### Reservation step-down

Before reserving the next stage, read Linux `MemAvailable` and reduce it by any
readable finite cgroup-v2 ancestor headroom (`memory.max - memory.current`).
Half of this available headroom is the reservation budget. Resident parents are
already included in current memory use; the other half remains available for
system activity, worker scratch and later allocations. If available memory
cannot be read, skip preallocation.

Compare the estimate against the actual power-of-two bucket capacities at the
87.5% load threshold. Repeatedly halve the estimate until the reservation fits.
If allocation still fails, release all partially reserved empty buckets, halve
again and retry; zero means ordinary on-demand growth. `try_reserve()` leaves
an existing linear-probe set intact on allocation failure. The original model
prediction remains unchanged for subsequent estimator calibration.

`GC_TRIANGLE_RESERVE_MIB=N` optionally lowers the preallocation budget; zero
disables preallocation. This is a reservation cap, **not a process memory limit**.
Subsequent insertion/growth still follows the existing policy and can exceed it.
Memory availability can change after the check; the policy is not an OOM guarantee.

Validation: `test-split-reservation` covers rounded capacities, step-down
boundaries, very large estimates and allocation failure under a child-only
address-space limit, checking that existing entries survive. Existing
linear-probe tests pass. Loop 9 with eight workers matches all graph bytes and
candidate counts with the automatic budget, a 1 MiB cap (six stages stepped down)
and zero reservation (ten stages). At V=13, the 1 MiB cap steps the 85,697-graph
estimate down to 10,712. Logs are in
`output/triangle_splits_L9_reservation_stepdown`.

### Latest loop-10 benchmark and loop-11 preparation

With adaptive estimates and reservation step-down enabled, three sequential
loop-10 runs with eight workers and `-O1` took 16.339, 22.900 and 22.614 seconds.
Median wall time was 22.614 seconds; median peak RSS was 285.0 MiB. The automatic
reservation budget was used, with no explicit cap. All 7,833,332 graph records
and 25,775,864 candidate counts match the saved reference in every run.
Compilation and comparison are excluded. Results and source snapshots are in
`output/triangle_splits_L10_current_benchmark`.

The current loop-11 OpenMP binary is built. Its 8,650 unrestricted seeds were
checked for dimensions, edge/loop counts, minimum degree three, triangle coverage,
connectedness and manifest counts. The prepared launcher is
`output/triangle_splits_L11_ready/run.py`; run it with Python 3 to record progress,
wall time and peak RSS. The output directory is fresh. This preparation did not
start loop 11; no completion or peak-memory claim is made for that loop.

### Producer experiments

Shared buckets are the selected implementation. The experimental producer and
ring-buffer code and build target have been removed; source snapshots and results
remain under `output/triangle_splits_L9_buffered_benchmark` and
`output/triangle_splits_L9_ring_benchmark`.

In the final paired loop-9 comparison with eight total threads, shared buckets
took 0.627 seconds median versus 1.307 seconds for the ring producer. Both
matched every graph file and candidate count. Loop 9 remains the quick benchmark
for further work.

### Optional OpenMP splitting

The serial executable remains the default. Build the separate OpenMP executable
with `make -f GraphGeneration/Makefile.pipeline triangle-splits-omp TRIANGLE_SPLIT_LOOP=10`.
Run it with `OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE ./build/triangle_splits_omp_L10
SEED_DIRECTORY NEW_OUTPUT_DIRECTORY` (on one shell line).

Workers expand different parents using the same cut-first splitting function.
Canonicalization stays local to each call. Workers insert into one shared
next-stage set, partitioned into 64 hash buckets with one mutex per bucket.
Only insertion and its candidate counter take the lock; canonicalization runs
outside it. Duplicate children are removed immediately, with no worker-local
sets or merge phase. A barrier precedes sorted output. Worker exceptions are
caught and rethrown on the main thread.

Loop-9 checks against saved serial output match every graph file byte-for-byte
and every candidate/output count for serial and 1, 2, and 4 OpenMP workers.
Measured times are 6.453, 6.541, 3.881, and 2.164 seconds respectively. Four
workers at loop 10 likewise match all 7,833,332 graphs and candidate counts,
taking 97.320 seconds versus the earlier serial 224.840 seconds. Peak RSS is
932,288 KiB (910.4 MiB). Results are under
`output/triangle_splits_L9_omp_workers{0,1,2,4}_run1` and
`output/triangle_splits_L10_omp_workers4_run1`, with sibling logs and metrics.
Compilation is excluded, and these runs use the same cut-vertex-inclusive seeds
from `output/triangle_comparison_L10_cut_splits_run1/seeds`.

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
`GraphStandardizer::standardize_no_sign` and deduplicated in the existing
`linear_probe_set`. Entries contain only the graph; an all-zero graph marks an
unused slot because every valid graph is simple. The stage uses 64 hash buckets
in both serial and OpenMP builds. Bucket selection uses high hash bits, leaving
the low bits available for probing. A resize temporarily duplicates only one
bucket rather than the entire stage. Before expansion, copy into a contiguous
vector and release each bucket immediately after copying it, then sort once by
canonical edges to preserve the previous file and parent order.
The splitting path uses neither `TransientGraph2` nor `BasisElement`/`LinComb`
operations.

The linear-probe replacement matches all 326,192 loop-9 graphs byte-for-byte
against the saved geng-validated run, including every stage's candidate counts.
Serial and OpenMP runs with 1 and 4 workers agree. Three serial runs with the
same `-O1` flags gave median wall times of 6.544 seconds for `std::set` and
5.925 seconds for `linear_probe_set` (9.5% less time). Peak RSS was approximately
18.8 MiB and 19.1 MiB respectively; this run shows no memory reduction. Sorting
briefly holds both the hash table and the compact vector. Compilation is excluded;
inputs are `output/triangle_comparison_L10_cut_splits_run1/seeds`. Measurements
are in `output/triangle_splits_L9_linear_probe.metrics.json`.


For each graph, scan all vertices and split every vertex of valence greater
than three. Enumerate subsets of incident edges to move to the new vertex,
leaving incidence zero on the old vertex to identify complementary splits.
Each side receives at least two original incidences and the connecting edge,
so both child vertices have valence at least three. Call `Graph::splitGraph`
to construct each plain child graph. Every split adds one vertex and one edge,
preserving loop order. There are no shared incidences or valence-sorting pass.

For cut-free parents, `SplitReductionRule.hpp` checks the proposed incidence
partition before constructing a child. Parent adjacency and degrees are cached
once. For a split set, update a temporary bitset adjacency view and reject if
another edge has a strictly larger sum of endpoint degrees and its contraction
gives a simple, cut-free graph. Endpoint common-neighbour and connectivity tests
are the same as in the earlier post-split filter. Ties survive; cut-vertex repairs
bypass the reduction filter.

Before enumerating split sets, skip vertex `v` entirely if the parent has an edge
`ab` disjoint from `v`, with no common neighbour, with
`degree(a) + degree(b) > degree(v) + 2`, and with the parent minus `a,b` connected.
The new edge's degree sum is always `degree(v) + 2`, regardless of the partition.
Contraction of this disjoint eligible edge commutes with splitting `v`, so it
rejects every partition of that vertex. Only surviving partitions call
`Graph::splitGraph`.

Why this preserves completeness: the alternate contraction has one fewer vertex
and edge, preserves loop order and minimum degree three, and is cut-free. Given
complete preceding stages, it is another available parent. Among all such
contraction edges of a child, a maximum-score edge survives the filter, so at
least one construction remains. Splits of cut-free parents remain cut-free;
therefore this filter does not change the internal cut-vertex repair stages.
No automorphism computation or additional canonicalization is needed.

The focused `test-split-reduction` target checks 28,255 splits against an
independent oracle that constructs contracted graphs and runs the existing cut
finder. It exhausts simple six-vertex parents of minimum degree three and adds
a separating-pair example and all 24 splits of K4,4, whose eight vertices can
all be skipped. Every tested partition agrees with the explicit-contraction
oracle, including vertices skipped as a whole. Loop 9 retains all 326,192 graphs byte-for-byte while
reducing candidates from 4,265,409 to 1,095,580 (74.3%). Two serial runs took
2.526 and 2.363 seconds versus the unfiltered linear-probe median of 5.925
seconds, with the same `-O1` flags and seeds. Four OpenMP workers took 0.793
seconds and matched serial output and candidate counts. Compilation is excluded.
Measurements are in `output/triangle_splits_L9_reduction.metrics.json`.

Exact checks through loop 10 preserve all 8,177,015 saved geng-validated graphs.
At loop 10 alone, candidates fall from 121,089,308 to 25,775,864 (78.7%),
with all 7,833,332 graph records byte-identical. The validation run took 63.080
seconds and peaked at 373.0 MiB RSS; some smaller-loop compilation overlapped,
so this is not an isolated timing comparison. Output and metrics are under
`output/triangle_splits_L10_reduction_run1`.


Moving the filter before child construction and skipping whole vertices improves
the shared-bucket loop-9 benchmark from 0.609 to 0.563 seconds median (7.5%).
Three alternating runs per variant use eight OpenMP workers, identical seeds
and `-O1` flags, with compilation excluded. All 326,192 graph records and
1,095,580 candidate counts are unchanged. Raw results and source snapshots are in
`output/triangle_splits_L9_pre_filter_benchmark`. Earlier loop-10 validation above
refers to the post-split implementation; this change was benchmarked on loop 9.

A leanness check compared whole-vertex skipping alone with whole-vertex plus
split-set filtering, keeping shared buckets and cut repairs in both. In three
alternating loop-9 runs with eight workers and `-O1`, vertices-only took 0.646
seconds median and canonicalized 1,521,273 candidates; the full filter took
0.564 seconds and canonicalized 1,095,580. Both produced identical graph files.
The per-partition check avoids 425,693 additional canonicalizations and is
retained. A split can break a triangle and make an edge eligible that was not
eligible in the parent, so the two checks are not equivalent. A concrete
same-vertex keep/reject example, sources and metrics are in
`output/triangle_splits_L9_vertices_only_benchmark`.

The current pre-construction filter with shared buckets was also timed at loop
10 with eight workers and the same `-O1` flags. Three sequential runs took
17.586, 19.038 and 19.322 seconds: median 19.038 seconds and median peak RSS
299.3 MiB (range 298.5–306.8 MiB). All 7,833,332 graph records and 25,775,864
candidate counts match the saved reference on every run. Compilation and
comparison are excluded. Results, logs and source snapshots are in
`output/triangle_splits_L10_pre_filter_benchmark`.

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
incoming splits that survive filtering and loaded seeds, before deduplication. Stage time includes
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
