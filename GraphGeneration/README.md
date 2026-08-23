# Graph representative generation

Build a generator specialized for one loop number:

```sh
make graph-stage-generator GRAPH_STAGE_GENERATOR_LOOP=7
```

Generate the vertex-irreducible stages from the unrooted one-edge support seed
`K2`:

```sh
OMP_NUM_THREADS=8 ./graph_stage_generator_7 generate output/L7
```

Transient files from the earlier rose, W3, and rerooting-capable generators
are not compatible. Use a fresh
output directory for a new run.
Commands that consume a transient file (`next`, `stats`, `validate`, and
`compare-transient`) reject the old payload kind.

For every vertex count `V`, generation writes two files:

- `loop_7_vertices_V_transient.gcg` contains canonical unrooted support
  families that have an active maximum-valence root allocation and can still be
  split. No root or marker hair is stored.
- `loop_7_vertices_V_admissible.gcg` contains canonical simple, hairless
  representatives together with their vertex-automorphism order and survival
  bits for the odd-edge and even-edge/odd-vertex sign conventions.

The one-stage commands are:

```sh
./graph_stage_generator_7 seed transient_V3.gcg admissible_V3.gcg
./graph_stage_generator_7 next transient_V3.gcg transient_V4.gcg admissible_V4.gcg
./graph_stage_generator_7 stats transient_V4.gcg admissible_V4.gcg
./graph_stage_generator_7 validate transient_V4.gcg admissible_V4.gcg
./graph_stage_generator_7 dimensions output/L7
```

`validate` checks that every transient support is unrootedly canonical and has
an active maximum-valence root allocation. It also
recomputes every final graph's canonical representative, automorphism order,
and sign-survival bits.

To compare a new final file with the simple subset of a legacy stage file:

```sh
./graph_stage_generator_7 compare-simple \
  output/L7/loop_7_vertices_8_admissible.gcg \
  reference/loop_7_vertices_8.gcg
```

For this generator, the relevant regression command filters the reference to
its vertex-irreducible subset:

```sh
./graph_stage_generator_7 compare-irreducible \
  output/L7/loop_7_vertices_8_admissible.gcg \
  reference/loop_7_vertices_8.gcg
```

`compare-transient` compares canonical unrooted-support sets produced by two
runs, which is useful for checking thread-count independence. It expects the
new support-transient format on both sides; it does not compare against old
one-hair frontier files.

## Generation rule

`unrooted_support_transient_graph<V,E>` wraps a loopless simple support graph.
Tadpoles are deleted and parallel bundles are collapsed to one support edge.
The fixed dimensions recover the total omitted surplus:

```text
surplus = E - number_of_support_edges
        = tadpoles + sum(bundle_multiplicity - 1).
```

A support record denotes all feasible allocations of that surplus rather than
remembering one tadpole count or one bundle-multiplicity vector. Only families
containing an active allocation enter the frontier. This merges states that
differ only by how tadpoles and excess parallel edges are distributed.

Generation starts at `V=2` from the unrooted one-edge support `K2`. At loop
number `L` it has `E=L+1` and therefore carries aggregate surplus `L`. The seed
permanent file is empty.

For one expansion, every vertex which can carry a feasible maximum-valence
exact allocation is considered as a temporary root. This need not be a vertex
of maximum support degree: hidden copies can raise a lower support-degree
vertex. The existing symbolic surplus rules enumerate which incident bundles
move and which are exhausted. The child then forgets that temporary choice and
is canonicalized as an unrooted support graph, so isomorphic root choices share
one frontier record.

Support leaves are allowed only next to the temporary root and require two
hidden copies of that bundle, raising their exact valence from one to three.
Likewise, a symbolic split may move no old support group when at least two
separated tadpole incidences stabilize the new exact vertex. Support bivalent
vertices require one hidden copy. Exact GC vertices are never univalent or
bivalent.

The augmentation is ordered by maximum support valence before canonicalization.
A split which creates a support leaf, or whose child has larger maximum support
valence than its parent, is emitted only when the parent already has a universal
support vertex, equivalently `maximum_support_valence == V-1`. This restriction
is tested by exact representative-set comparison through loop number eight.

The unrooted standardizer follows the `standardize4` refinement scheme. Its
first colors are literal support valences: it sorts and groups them before any
hashing. Maximum-support vertices therefore form the final initial cell. The
standardizer returns both its valence and its size. In memory the registry is
partitioned by `(V,E,maximum_support_valence,tied_maximum_vertices)` and then
sharded within each partition.

The retained root side receives strictly more exact old half-edges than the
diverging side for every continuing split. Equality is terminal. In support
degree the two sides may tie only when a selected bundle remains on both sides
(a virtual double edge). Support-leaf splits, which select no old root bundle,
are allowed only after the support has reached a universal vertex; the
temporary root itself need not be that vertex.

An admissible image must have zero surplus, so its support is the exact simple
graph. It is finalized independently of the raw parent-root label and
canonicalized without a distinguished root before insertion into the
permanent set. This second canonicalization is mandatory: forgetting the root
does not preserve the rooted canonical labeling.
Terminal simple images enter only the permanent set; active simple images
enter both the permanent set and the next transient frontier.

The transient registry is keyed by canonical unrooted supports.
There is no construction spanning tree, preferred raw witness, stored root
hair, tadpole count, or parallel-bundle multiplicity vector.

This unrooted maximum-choice construction remains experimental. Compare its
permanent representative sets against the legacy vertex-irreducible sector
before treating new dimensions as complete.

For loop number `L`, a complete generation would end at
`V = 2(L - 1)`. At that point
`E = L + V - 1` and minimum valence three force every graph to be trivalent;
the final transient file is therefore empty while the permanent file is kept.

## Dimensions and progress

After closing each file, the generator prints and flushes its kind, loop number,
vertex count, edge count, number of records, and cumulative stage time.
Permanent-file lines also report both signed survivor counts.

Each loop directory receives `generation_dimensions.tsv` with columns:

```text
loop  vertices  edges  unsigned_classes  odd_edge_classes  even_edge_odd_vertex_classes
```

Transient counts are intentionally excluded from this public dimension table;
they are only implementation diagnostics.

`stats` reports support-family and admissible counts. It no longer reports a
transient valence range because one support family can represent several exact
valence vectors.

To generate and combine several loop orders:

```sh
OMP_NUM_THREADS=8 make graph-dimension-table \
  MIN_TADPOLES=3 MAX_TADPOLES=11
```

The default maximum remains 5 so an unexpectedly large computation is not
started accidentally.

## Mapped file format

Version-2 mapped files begin with a 64-byte `MappedGraphFileHeader` containing
the graph dimensions, record count, record stride, and payload kind.

Ordered unrooted support records with repairable leaves use payload kind 12.
The unordered repairable-leaf experiment used kind 11. The
leaf-free unrooted experiment used kind 10. K3-seeded immutable-root
ordered-split records use kind 9. The prior W3-seeded immutable-root records
use kind 8, rerooting-capable records use kind 7,
minimum-support-degree-three records use kind 6, and the older
K3/minimum-support-degree-two records use kind 5; all are rejected. The
`V(V-1)/2` possible support
edges are packed consecutively, one bit per edge in lexicographic vertex-pair
order, so the record stride is

```text
ceil((V choose 2) / 8) bytes.
```

The expanded edge count `E` comes from the header, so neither the omitted
surplus nor any of its allocations is written. Vertex labels carry no root.
Every stored support record is splittable, hence the header's splittable count
equals its total record count.

Final records contain the `2E` endpoint bytes of a canonical hairless graph,
followed by a `uint64_t` vertex-automorphism order and one byte of sign-survival
flags. Metadata is copied with `memcpy`, so records need no alignment padding.

Record order is deliberately unspecified because the parallel registries are
unordered. Consumers must build basis indices from the file they actually
read; indices are not promised to remain stable across separate generation
runs.

Legacy version-0 and version-1 graph-only files remain readable as hairless
reference inputs. Old version-2 one-hair and unrestricted-support frontier
files use different payload kinds and are rejected as transient inputs,
preventing incompatible generation sectors from being mixed silently.

## Focused tests

```sh
make run-transient-graph-test
make run-support-transient-test
make run-mapped-support-transient-test
make run-final-canonicalization-test
make run-sparse-rank-test
```

## Performance comparison with geng

The generator and geng benchmarks are deliberately separate runtime scripts;
the graph generator has no nauty headers or libraries in its build dependency
graph. Build the specialized generator first, then run the two measurements
into the same fresh benchmark directory and compare their counts and timings:

```sh
make graph-stage-generator GRAPH_STAGE_GENERATOR_LOOP=8
tools/benchmark_gc_generation.sh 8 4 /tmp/gc-geng-L8
tools/benchmark_geng_generation.sh 8 /tmp/gc-geng-L8
tools/compare_generation_benchmarks.sh /tmp/gc-geng-L8
```

The geng side discovers `nauty-geng` or `geng` at runtime; `GENG=/path/to/geng`
overrides discovery. It runs `geng -C -d3 -u` at the single edge count
`E=V+L-1`, matching the biconnected minimum-degree-three permanent sector.
The comparison refuses to report a timing ratio unless every `(V,E)` count
agrees. Compilation time is excluded. GC peak memory is measured for the whole
multi-stage run; geng is launched once per vertex count and its reported total
is the sum of those stage wall times, with peak memory the largest stage.

The support test exercises quotient generation against exact multiplicity
generation from W3. The mapped-support test
covers compact round trips and rejection of incompatible headers and record
layouts. The final canonicalization test cross-checks the stored survival logic
against the two existing signed standardizers.

For a comparison closer to producing the actual graph-complex basis, build the
stream finalizer and run geng without `-u`:

```sh
make geng-finalizer GRAPH_STAGE_GENERATOR_LOOP=8
tools/benchmark_geng_finalized_generation.sh 8 /tmp/gc-geng-L8
tools/compare_finalized_generation_benchmarks.sh /tmp/gc-geng-L8
```

This streams every graph6 record from geng through the same final
canonicalization used by the GC generator, computes the vertex-automorphism
order and both sign-survival flags, and writes compatible final `.gcg` stage
files.  It remains a one-worker pipeline and retains no stage-wide graph set.
