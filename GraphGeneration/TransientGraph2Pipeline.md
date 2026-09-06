# New graph generation pipeline

These are the current generation rules for `transient_graph2` and the
single-threaded `GraphGenerationPipeline`. Completeness has not yet been proved.

## Build and run

The new pipeline has its own build file; the old generator is unchanged.
From the repository root:

```sh
make -f GraphGeneration/Makefile.pipeline PIPELINE_MAX_LOOP=3 PIPELINE_MAX_VERTICES=6
./build/transient_pipeline_L3_V6 output/transient2_L3
make -f GraphGeneration/Makefile.pipeline test
```

### Counts and performance against geng

Build the comparison's generator through loop 6 and 10 vertices:

```sh
make -f GraphGeneration/Makefile.pipeline comparison \
  PIPELINE_MAX_LOOP=6
```

Run TransientGraph2 first, then geng separately, and compare counts:

```sh
make -f GraphGeneration/Makefile.pipeline run-comparison \
  PIPELINE_MAX_LOOP=6
```

For a three-way comparison using the existing hairless implementation:

```sh
make -f GraphGeneration/Makefile.pipeline \
  run-comparison-three-way \
  PIPELINE_MAX_LOOP=6 \
  COMPARISON_DIR=output/three_way_L6_run1
```

`comparison-three-way` builds without running. It reuses
`tools/generate_graph_stages.cpp` and its current unrooted support rules,
compiled into `build/hairless_pipeline_L3` through the selected loop bound.
This is the existing implementation as it stands, not a reconstructed historical
revision or an assumption that it is correct. Both generators use `PIPELINE_FLAGS`;
the old generator additionally requires OpenMP.

To run only hairless versus geng through loop 8:

```sh
make -f GraphGeneration/Makefile.pipeline \
  run-comparison-hairless \
  PIPELINE_MAX_LOOP=8 \
  COMPARISON_DIR=output/hairless_vs_geng_L8_run1
```

This builds and runs only the hairless generators and geng. The final count
table lists loop order, vertex and edge counts, each backend's final graph
count, and match status. Totals and timings are printed separately. The tables
are also saved in `reports/counts.txt` and `reports/performance.txt`.

The runner executes TransientGraph2 once through the maximum loop order, then
the old hairless generator separately for each exact loop order, then geng.
Hairless uses one thread by default; set `HAIRLESS_THREADS` to change it.
TransientGraph2 currently uses minimum split valence 2 with only the triangle seed. The former
`GC_DISABLE_BIVALENT_SPLITS` experiment switch is retired.

The extra `hairless` and `hairless_status` columns in `counts.tsv` compare unsigned
2-connected minimum-degree-three counts from the old generator's
`generation_dimensions.tsv` against geng. Its sign-survival counts are not used.
Files and live logs are under `hairless/L<loop>` and `hairless/L<loop>.log`.
Failures and missing rows are reported explicitly, never converted to zero;
an old-generator failure does not prevent comparing the other backends.

`performance.tsv` includes each exact-loop hairless run and their total.
Failed totals are marked incomplete. The old pipeline also computes automorphism
orders and sign metadata, so these are end-to-end workflow timings, not timings
of identical work. Counts alone do not prove equality of graph sets.

The comparison target automatically uses `2 * (PIPELINE_MAX_LOOP - 1)` vertices.
The shell runner requires nauty's `geng` (`pacman` package `nauty`), and uses
Bash timing and standard shell utilities. It uses no Python.

The default output is `output/transient2_counts_L6`. It must not already exist;
set `COMPARISON_DIR=output/another_run` for subsequent runs. The directory contains:

- `transient/`: pipeline stage files.
- `transient.log`: full generator output; the terminal shows labelled progress only.
- `geng/`: reference graph6 files, stderr logs, and per-stage timing files.
- `counts.tsv`: both counts and match/mismatch status for each feasible `(L,V,E)`.
- `performance.tsv`: wall, user CPU, and system CPU seconds for the pipeline,
  each geng invocation, and their summed reference time.
- `reports/counts.txt`: a count-only summary for each backend, with match,
  mismatch, and unavailable row counts against geng.
- `reports/performance.txt`: a separate timing-only summary for each backend.
- `reports/counts/{transient2,hairless,geng}.tsv`: individual backend count
  tables by loop, vertices, and edges (hairless only in three-way runs).
- `reports/performance/{transient2,hairless,geng}.tsv`: individual backend
  timing tables. Failed runs and incomplete totals are labelled.

To produce these reports from an existing run without generating anything:

```sh
bash tools/compare_transient2_counts.sh \
  --report output/three_way_L7_run1
```

The count summary excludes intermediate transient graphs. Its match counts
refer to partition rows, not individual graphs; graph totals can be partial
if results are unavailable. Wall time means elapsed time, while user and system
times measure CPU usage. Full diagnostic logs remain separate from these reports.

The generator now prints each saved partition immediately. While expanding a
stage, it reports completed parents about every five seconds, checked between
parents; one expensive parent can delay the next report. Stage timings include
writing that stage and expanding it into the next. An already running older
binary will retain its previous logging behavior.

The comparison covers vertex-2-connected simple graphs of minimum valence three,
through the selected loop number. Pipeline counts apply the valence filter and
require connectivity after deleting each vertex while writing the existing
stage files; the reference uses `geng -C -d3` at the
same vertex and edge counts. Matching counts alone do not prove that the graph
sets agree, and this does not check odd/even GC signs. A missing count, failed
generator, or count mismatch makes the runner exit nonzero.

Transient files still retain intermediate graphs with cut vertices; a later
split can remove the cut vertex. The filter applies to reported GC counts and
to the separate transient-to-GC converter. The runner checks the generator's
`COUNT_SCOPE biconnected_min_degree_3` marker to reject older count semantics.
Use a fresh output directory when comparing this scope with geng. Here
"3-valent" means minimum valence three, not exactly cubic.

Compilation is excluded from the timing. Both workflows include file output;
the pipeline additionally constructs and writes intermediate low-valence
graphs. These are end-to-end workflow timings, not identical workloads.

The equivalent C++ interface is:

```cpp
GraphGeneration::GraphGenerationPipeline<3, 6> pipeline;
auto summaries = pipeline.run("output/transient2_L3");
```

The loop and vertex bounds are compile-time parameters because graph dimensions
are part of the C++ type. A pipeline instance runs once and refuses to overwrite
existing stage files. The triangle seed requires a loop bound of at least one.

## Current generation rules

- Store connected simple graphs, partitioned by vertex count and loop order
  `L = E - V + 1`.
- Seed only with the triangle `K3`. No star seeds are inserted.
- For every standardized parent, `BivalentSplitRule.hpp` scans the low-valence
  prefix, skips leaves, and selects only vertices adjacent to every bivalent
  vertex. Bivalent vertices themselves cannot be selected. With none present,
  every vertex is eligible. The triangle seed is the sole exception: split
  just vertex 0, since its three vertices are equivalent.
- Each incident edge goes to either split vertex or both. Connect the two
  split vertices. Both must have valence at least **2**, including that edge.
  At least one split vertex must also reach the parent’s second-largest
  valence, `valences[V - 2]`, including the connecting edge. The same threshold
  applies to every vertex selected in that parent. Both split vertices must
  also stay at or below `valences.back() + (valences.back() == V - 1)`: the
  parent maximum plus one if any vertex is universal, otherwise the parent
  maximum. This cap applies to the two split vertices. Enumeration prunes a
  branch as soon as either side exceeds it.
  Starting from K3, these splits cannot produce leaves. The separate
  extended-edge operation is currently disabled in the pipeline.
- Each shared incidence adds one edge and one loop. The pipeline owns the
  partition budget `max_surplus_edges = MaxLoopNumber - L`.
- The split API is
  `split(vertex, preserve_valence, min_valence, max_valence, min_loop_number, max_loop_number, collector)`.
  When the selected vertex has bivalent neighbours, the pipeline calls
  `split_with_bivalent_neighbours` with the same arguments. This separate
  entry point checks completed incidence assignments before constructing a
  child: either no bivalent vertices remain, or the split pair includes a
  bivalent vertex and one non-bivalent vertex in the child neighbours all
  bivalent vertices. Thus a split producing two vertices of valence at least
  three must eliminate all bivalent valencies from the child. Assignments failing this test never reach the collector.
  Ordinary `split` does not apply this test.
  The pipeline passes minimum valence 2, minimum loop order 0, and maximum
  loop order `L + max_surplus_edges`. Both loop bounds are inclusive.
  A higher minimum loop order requires at least
  `max(0, min_loop_number - L)` shared incidences; branches that cannot reach
  that minimum are pruned by the split operation.
- Standardize every child without signs and deduplicate it in its partition.
  **There is no child eligibility filter.** The common-neighbour rule only
  selects parent vertices for splitting; the specialized split additionally
  rejects assignments without an eligible child root before child construction.
- Standardization orders valencies increasingly. Vertex selection does not
  require maximum valence, and there is no surplus-allocation selection rule.
- Exchanging the two split vertices is suppressed during split enumeration;
  canonicalization resolves the remaining isomorphisms.
- Process vertex stages through `MaxVertices`; do not expand the last stage.

The minimum-valence-three and vertex-2-connectivity checks apply only to the
reported comparison counts and the separate GC converter. Every emitted transient is standardized and stored, including graphs with
leaves, bivalent vertices, or cut vertices.

## Separate bivalent-creation operation

`TransientGraph2` provides
`create_bivalent_vertices(vertex, preserve_valence, max_valence, min_loop_number, max_loop_number, collector)`.
For each edge `vertex--v`, it emits one graph retaining that edge and adding
`vertex--b--v`, where `b` is a new bivalent vertex. It adds one vertex, two
edges, and one loop. The child loop number must lie within the inclusive
minimum and maximum; otherwise nothing is emitted. Before constructing any
child, require `preserve_valence <= degree(vertex) + 1 <= max_valence`. Only
the chosen vertex is checked; the neighbour and new vertex are not checked
against these valence bounds. Isomorphic results from
different incident edges are left to the caller to standardize and deduplicate.
This operation remains implemented, but the pipeline currently does not call
it. The current benchmark uses minimum split valence two and the triangle
seed alone; it retains the existing vertex-selection, preserve-valence,
maximum-valence and specialized bivalent-assignment rules.
Ordinary pipeline splits produce vertices of minimum valence two;
the generic split primitive still accepts a caller-supplied minimum for testing
and reuse. Stored seed vertices may have lower valence.

## Storage and progress

The pipeline is single-threaded. It collects children in typed buffers,
standardizes and deduplicates them after each parent, writes completed stages,
and releases previous-stage storage after expansion.

Files are named `transient_V<V>_E<E>.gcg`, including empty partitions within
bounds. Each has the existing 64-byte version-2 header with payload identifier
13, followed by `2 * E` endpoint bytes per graph. Both header count fields
hold the stored record count. Files are written through a temporary file and
renamed; an existing output file or temporary file is never overwritten.

Logs report all transient records separately from the GC comparison counts,
with stage timings and parent progress approximately every five seconds.

## Validation

The integration test covers bivalent root selection, the single triangle-seed
split, and rejection of the C4 assignment before constructing a child,
the absence of star seeds at larger loop budgets, canonicalization, deduplication, file
round-trips, repeatable output, and invalid files. The split unit test compares
against exhaustive incidence assignments across minimum and maximum loop bounds.

Earlier matching geng counts through loop 8 used restricted vertex selection;
they do not validate the current bivalent selection and preserve-valence policy. Use a fresh comparison directory.
