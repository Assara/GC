# Triangle-seed graph generation

At loop number L, the pipeline generates connected simple graphs of minimum
valence three. A seed has every edge covered by a triangle. The unrestricted
seed generator starts with K3 and adds triangles with zero, one, or two new
vertices, identifying graphs by canonical edges alone. Its completion budget
prunes graphs that cannot reach the requested loop cap with minimum valence three.

```sh
make -f GraphGeneration/Makefile.pipeline triangle-seeds TRIANGLE_MAX_LOOP=9
./build/triangle_seeds_unrestricted_L9 output/triangle_seeds_L9_new
```

Seed generation writes a counts.tsv manifest and seeds_L*_V*.g6 files. Seeds
include cut vertices: the splitting stage resolves cuts before expanding other
vertices and emits only cut-free graphs. Keep these intermediate seeds when
splitting; filtering them out loses valid output graphs.

## Splitting and filtering

Splitting degree d moves k incidences to a new vertex and adds the connecting
edge, giving degrees k+1 and d-k+1. Both sides must have degree at least three.
Incidence zero remains on the old vertex to eliminate complementary duplicates.
Loop number stays fixed; vertex stages end at V=2(L-1).

For a cut-free parent, SplitReductionRule skips a whole vertex when a disjoint
triangle-free edge has a larger endpoint-degree sum than the proposed new edge
and deleting its endpoints leaves the parent connected. A second check applies
the same preference to virtual child adjacency before constructing the child.
Equal scores remain eligible. Cut repair bypasses this reduction and uses
CutVertexSplitRule. The tests compare the reduction against actual contractions.

Both splitters require a fresh output directory and save graph6 files and
counts.tsv. The shared-bucket splitter keeps a complete vertex stage in RAM.
The valence-group splitter uses disk-backed parent groups to reduce memory.

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

## Shared-bucket implementation

```sh
make -f GraphGeneration/Makefile.pipeline triangle-splits-omp TRIANGLE_SPLIT_LOOP=9
OMP_NUM_THREADS=8 ./build/triangle_splits_omp_L9 \
  output/triangle_seeds_L9_new output/triangle_splits_L9_new
```

For serial execution, build `triangle-splits` and use `build/triangle_splits_L9`.

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

## Validation

```sh
make -f GraphGeneration/Makefile.pipeline test
make -f GraphGeneration/Makefile.pipeline build/valence_split_pipeline_test_L9
./build/valence_split_pipeline_test_L9 REFERENCE_DIRECTORY VALENCE_OUTPUT_DIRECTORY
```

The comparison test verifies exact canonical graph identities, absence of
duplicates, sorted valences, relabelling invariance, and split routing.

For independent comparison with nauty geng and labelg (installed separately):

```sh
make -f GraphGeneration/Makefile.pipeline run-triangle-comparison TRIANGLE_MAX_LOOP=6
```

This generates seeds and all split stages through the cap and compares each
stage to `geng -C -d3`. The runner supports caps 3 through 10. Missing, extra,
and duplicate isomorphism classes cause failure. Use fresh output directories.

## Completed large run

The L11 shared-bucket run saved 216,162,805 graphs in 533.55 seconds with eight
workers, peaking at 5.93 GiB RSS. Its graph6 files occupy about 4.3 GiB under
`output/triangle_splits_L11_ready/graphs`. The final cubic stage has 497,818
graphs. Completion was successful; unlike L10, this run has not had an independent
full graph-identity comparison. L12 has not been run through splitting.
