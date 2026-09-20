# Enumerated bases and contraction matrices

The first homology stage loads generated graph6 files for fixed loop L and
vertex counts V-1, V, V+1, removes graphs that vanish under the chosen orientation
signs, and enumerates canonical representatives in lexicographic half-edge order.
All three bases and the two contraction matrices stay in RAM in a
`GraphHomology::ContractionWindow<L,V,Parity>`.

```sh
make -f GraphHomology/Makefile HOMOLOGY_LOOP=6 HOMOLOGY_VERTICES=8
./build/contraction_matrices_L6_V8 GRAPH_DIRECTORY even
```

Loop and vertex counts specialize the executable; parity is selected at runtime.
`even` means odd edges and even vertices (the wheel convention, c=0,d=1).
`odd` means even edges, odd vertices and antisymmetric edge directions (c=d=1).
Coefficients are interpreted in the current `fieldType`, presently F_32783.

`down` has rows indexed by the V-1 basis and columns by the V basis. `up` has
rows indexed by V and columns by V+1. They use the existing unscaled contraction
differential, including its cut-vertex rule. Parallel-edge targets are zero for both parities:
parent edges in triangles are skipped before constructing or canonicalizing the
child. Contributions are signed units;
canonical equal targets are combined and zeros removed before basis lookup.

Finished matrices and stored transposes use `VectorSpace::OwnedArray<T>`:
a move-only `unique_ptr<T[]>` wrapper that owns its length and provides
`size()`, `data()`, indexing and iteration. There is no resizing or spare
capacity. Nonzero and column counts are derived from the arrays.
Contraction construction uses a vector-backed builder while counts are unknown;
`finish()` copies into exact-size arrays, releasing each builder array after
copying it. This temporarily requires both copies of the array being converted.
Transposes allocate fixed arrays directly from the known dimensions.

Sparse storage consists of separate row-index (uint32_t), coefficient (int8_t),
and column-offset (size_t) arrays. This costs 5 bytes per nonzero plus one offset
per column and a final offset. Contributions accumulate in wider integers before
storage. With at most 127 source edges, signed-byte coefficients suffice; no
runtime coefficient-bound check is performed. This encoding applies to the
original differential, not to coefficients created by later field elimination.

Bases are stored once as packed edge endpoints. Matrix construction uses binary
search in the sorted target basis, avoiding a second graph hash table. IDs do not
depend on the order of records in input graph6 files. Duplicate input classes and
missing nonzero targets are errors. Input files must contain the generated
cut-free, simple graphs of minimum degree three. Mathematically impossible simple
graph stages at either end of the chosen window are treated as empty.

The diagnostic executable prints basis dimensions, nonzero counts, allocated
array bytes (not process RSS), and construction time. It checks the entire product
`down * up` over the integers without allocating a product matrix. The library
also supplies matrix-vector multiplication over the selected field. Matrices are
not written to disk. Optional block-Wiedemann rank computation gives homology
dimension estimates without extracting class representatives.

## Simple-graph convention

Both parities use the simple-graph quotient: contractions producing parallel
edges are set to zero. This is explicit for odd parity; in even parity such
graphs already vanish by the edge-orientation signs. The two-vertex, three-edge
multigraph is outside this basis and is not added as an exceptional class here.
The relationship to the full multigraph complex is not certified by this tool.

## Checks

```sh
make -f GraphHomology/Makefile test GRAPH_DIRECTORY=PATH_TO_L6_GRAPHS
```

The tests check signed-byte storage, cancellation, multiplication over the current
field, detection of a nonzero differential product, and every window in both parities
V=5 through V=10 at L6 against the existing graph contraction routine followed by an independent check
for parallel edges in the actual children. The tests verify that nonzero
parallel-edge terms are encountered and discarded. Both
matrices in every window satisfy d squared equals zero over the integers.

## Higher-loop construction checks

Single-threaded default -O1 builds were checked at (L,V) = (7,10), (8,11), and
(9,13), choosing the peak generated graph stage at each loop. Both parities pass
the full integer check of the product of consecutive matrices. At L9,V13:

| Parity | Dimensions V12 / V13 / V14 | NNZ down / up | Build seconds | Peak RSS |
|---|---|---|---|---|
| Even | 65,390 / 68,987 / 44,434 | 874,482 / 642,540 | 3.835 | 26.43 MiB |
| Odd | 66,525 / 68,725 / 43,423 | 880,062 / 633,798 | 3.744 | 25.94 MiB |

These are construction checks of selected windows, not homology calculations or
all-degree sweeps. Compilation is excluded; build seconds exclude the product
check. Peak RSS includes bases, matrices, allocation overhead and verification.
Full results are in `output/contraction_matrices_L7_L9_9a95xvkb`.

## Homology dimensions over the current field

```sh
make -f GraphHomology/Makefile HOMOLOGY_LOOP=7 HOMOLOGY_VERTICES=8
./build/contraction_matrices_L7_V8 GRAPH_DIRECTORY even --rank 8
./build/contraction_matrices_L7_V8 GRAPH_DIRECTORY odd --rank 8
```

The optional block width defaults to 8. After checking the chain relation, the
tool computes both ranks using three independently preconditioned block-Wiedemann
trials and reports `dim C_V - rank(down) - rank(up)`. The ranks and homology
dimension are explicitly labelled probabilistic estimates. Small L6 tests compare
all ranks to deterministic sparse elimination. Class representatives are not
computed. See [BlockWiedemann.md](../VectorSpace/BlockWiedemann.md) for the algorithm,
API, validation, and current performance limitations.

At L7,V8, the default -O1 run gives:

| Parity | dim C8 | rank down | rank up | Estimated dim H | Rank seconds |
|---|---:|---:|---:|---:|---:|
| Even | 75 | 10 | 64 | 1 | 0.0073 |
| Odd | 83 | 15 | 68 | 0 | 0.0073 |

All three rank trials agreed for both matrices. These timings exclude construction
and compilation. Logs are under `output/block_wiedemann_validation`.

## Capacity and literature comparison

The complete loop 3–8 sweep over F_32783 matches all 58 feasible (loop, vertex,
parity) entries of Figure 1 in Brun–Willwacher,
[Graph homology computations](https://nyjm.albany.edu/j/2024/30-5v.pdf).
Their simple, cut-free graph convention agrees with this pipeline. A vertex
count V at loop L has homological degree V-1-(n-1)L in G_n; choose n=2
for the even convention and n=3 for the odd convention.

All integer chain checks passed, all three rank trials agreed, and ranks agreed
across 46 overlapping windows. Loop 8 required 38.58 seconds of rank work
for both parities combined, with block width 8 and the default single-threaded
-O1 build. Adjacent windows currently recompute their shared differential rank.
These are probabilistic finite-field results, not rational rank certificates.

Raw logs, timings, reference entries, comparison code and loop-9 capacity
probes are saved in `output/homology_capacity/README.md`.

## Representatives from the transposed natural composition

```sh
make -f GraphHomology/Makefile HOMOLOGY_LOOP=8 HOMOLOGY_VERTICES=12
./build/contraction_matrices_L8_V12 GRAPH_DIRECTORY even --nullspace OUTPUT_DIRECTORY 8
```

`--nullspace` now uses `(S_down C_down + C_up S_up)^T` directly for Krylov
multiplication. Natural adjoints are constructed by exact automorphism scaling.
There is no additional Gram product. Integer accumulation uses `GraphAccumulator`
(default int32, switchable to int64), with reduction after both branches. A random
nonzero left diagonal is applied after reduction to precondition the recurrence;
it preserves the right kernel without changing the integer accumulation bound.
Narrow reconstruction blocks reuse the Krylov workspace.

Before saving, each vector z is converted to original contraction coordinates:
`x[G] = z[G] / |Aut(G)|`, and its first nonzero coefficient is normalized to one.
Both `C_down x = 0` and `S_up x = 0` are checked explicitly. A failure aborts
without exporting representatives. All automorphism orders must be units in the
field; this is checked before solving. There is no automatic field-retry logic.

The output directory must be new and its parent must exist. It contains:

- `basis.tsv`: zero-based graph IDs followed by ordered, directed edges specifying
  the stored orientation.
- `vectors.tsv`: vector ID, graph ID and nonzero coefficient in `0..p-1`, in
  original contraction coordinates.
- `automorphisms.tsv`: graph ID and integer automorphism order.
- `metadata.txt`: operator, conversion, verified residuals, field and solver settings.

The extracted vectors are independent; completeness depends on the probabilistic
rank estimate. Production does not independently certify their independence
modulo boundaries. Small L6 tests do check both the homology dimension and
independence modulo the image by exact elimination, in both parities.

`--legacy-nullspace OUTPUT_DIRECTORY [BLOCK_SIZE [SEED]]` retains the previous
unweighted `[down; up^T]` solver and its original output convention. Existing
compiled binaries and saved runs are not changed by this source update.

`--natural-adjoint [BLOCK_SIZE [SEED]]` constructs the natural adjoints by exact
automorphism scaling, checks the integer accumulator bound and evaluates one
block of the transposed natural composition. See
[NaturalCompositionBounds.md](NaturalCompositionBounds.md#implementation) for
the formula, 32/64-bit switch and validation. This construction mode does not
run nullspace extraction.
