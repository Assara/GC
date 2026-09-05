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

The equivalent C++ interface is:

```cpp
GraphGeneration::GraphGenerationPipeline<3, 6> pipeline;
auto summaries = pipeline.run("output/transient2_L3");
```

The loop and vertex bounds are compile-time parameters because graph dimensions
are part of the C++ type. A pipeline instance runs once and refuses to overwrite
existing stage files. The triangle seed requires a loop bound of at least one.

## Representation and stages

Transient graphs are connected and simple; ordinary 1-valent vertices are
allowed, with no external half-edge slots. Valence counts ordinary edge
incidences. The loop number is `L = E - V + 1`.

For a target maximum loop number `N`, a stage contains standardized,
deduplicated graphs with exactly `V` vertices and loop number at most `N`.
Split the stage's graphs, standardize and deduplicate their children in the
next collector, then write that completed `V + 1` stage to a file. Use the
completed collector as the input to the following stage.

Children with different edge counts have different concrete graph types.
The pipeline owns temporary, typed collector buffers and a `linear_probe_set`
for each `(V,E)` pair. After splitting each parent, it drains the buffers:
discard children with no eligible split vertex, standardize each remaining
child, insert it into the matching set, then clear the buffers. Candidate
counts include the children rejected by this filter.
Files are written directly from the sets, skipping empty hash-table slots.
Previous-stage storage is released after expansion. There are no threads,
locks, sharding, or separate admissible/final outputs in this pipeline.

## Transient file format

Each stage has files named `transient_V<V>_E<E>.gcg`, including empty partitions
within the configured bounds. Only transient graphs are written.

`MappedTransientGraph2File.hpp` adapts the existing mmap reader/writer without
changing the old implementation. Each file contains:

- The existing 64-byte version-2 header with `(V,E)`, record count, record size,
  and new payload identifier **13**. Both count fields equal the record count;
  flags are zero. This identifies stage records, without an active/terminal
  classification.
- Exactly `2 * E` bytes per graph: the standardized endpoint array in edge
  order, one byte per vertex label. No C++ object padding, pointers, hash-table
  capacity, or empty slots are serialized. The header retains the existing
  native integer layout; the endpoint payload itself is byte-based.

The writer maps a temporary file, fills it, closes the mapping, and renames it
to the stage filename. The new reader validates the payload identifier in
addition to the existing dimension, record-size, version, and file-size checks.

## Valency array and split parameters

The pipeline computes the parent's valency array. Standardization is expected
to order vertices by nondecreasing valence:

```text
valency_array[0] <= ... <= valency_array[V - 1]
```

The new standardizer canonizes first, then stably orders those canonical vertex
labels by valence and sorts the edges again. This enforces the ordering without
altering the old standardizer. The last vertex has index `V - 1`; the entry
before it (`last - 1`) is `valency_array[V - 2]`.

For each parent:

- Set `max_loop_number = N`.
- Set `preserve_valence = valency_array[V - 2]`, the parent's
  second-largest valence. For the five-vertex wheel this is 3, not 4.
- Set `min_valence = 1` for every parent.

Only split maximum-valence vertices adjacent to every other vertex of valence
1 or 2. These vertices form one combined group: separate common neighbours
for leaves and bivalent vertices do not suffice. The splitting vertex itself
is excluded from this adjacency requirement, allowing the triangle seed to
split. With no vertices of valence 1 or 2, every maximum-valence vertex is
eligible. One edge scan counts neighbours in the combined group.

The same eligibility check filters children before standardization and
deduplication. A child with no eligible vertex is discarded. Raw child labels
are not yet ordered, so the check computes the maximum valence explicitly.
There is no count cap or universal-preservation override. The separate GC
conversion requires minimum valence 3, discarding both leaves and bivalent
vertices. These changes allow transient leaves; the split operation itself
already supports `min_valence = 1`.

For the first implementation, the provisional defaults are: start from the
triangle `K3`, split every eligible maximum-valence vertex, and stop at the explicit
maximum vertex count.
These defaults still need to be justified alongside the generation rules.

## Rules inside `split`

```cpp
split(vertex, preserve_valence, min_valence, max_loop_number, collector);
```

- Each incident edge can attach to either split vertex or to both. Add an
  edge between the split vertices.
- Both split vertices must have valence at least `min_valence`; at least one
  must have valence at least `preserve_valence`.
- A split vertex may exceed the parent's splitting-vertex valence only if
  that parent vertex is universal. Universality is determined from its
  incident-edge count; it is not a parameter.
- Sharing `K` incident edges produces `V + 1` vertices, `E + 1 + K` edges,
  and loop number `L + K`. Do not exceed `max_loop_number`. No minimum loop
  bound is currently implemented.
- Fix the first exclusively assigned incidence on the original vertex to
  avoid mirror assignments. Shared incidences do not distinguish the sides;
  the all-shared assignment is emitted once. Other isomorphisms are resolved
  by standardization and collector deduplication.

There is no numerical cap on 1- or 2-valent vertices. The common-neighbour condition
is applied by the pipeline before standardization and when choosing a parent
vertex, not inside `split`.

## Proof to follow

We still need to prove that suitable seeds, vertex choices, and the stated
restrictions generate every graph in the intended target class up to loop
number `N`. Small-case comparisons can support that work but do not replace
a completeness proof.

The small integration test uses loop bound 3 and vertex bound 6. It retains
the triangle with a leaf and checks the combined common-neighbour condition
on every stored graph. Tests also check deduplication, canonical
valence ordering under vertex relabeling, file round-trips, repeatable bytes,
empty partitions, and rejection of incompatible/truncated files. Transient
stages are intermediate graph collections, not the final GC bases.

With the combined 1-/2-valent rule, the loop-7 check recovers all four formerly
missing graphs. It produces 1108 connected GC graphs, with 588 odd and 579
even survivors, matching the full connected reference. Graphs with cut
vertices are now also produced. See
[`output/transient2_L7_min1/REPORT.md`](../output/transient2_L7_min1/REPORT.md).
