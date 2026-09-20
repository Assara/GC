# Block Wiedemann rank, nullspace and solve

`block_wiedemann_solver<K>` computes rank/nullity over a finite field and solves
rectangular systems. Native signed-byte CSC matrices are borrowed directly:
the original is retained in signed-byte form and one signed-byte sparse transpose
is allocated. Both products gather into independent output rows using OpenMP
(default eight threads), following the scalar solver design. The existing
compressed matrix and legacy LIL entry points are supported (LIL owns one
compressed copy plus its transpose). Borrowed matrices must outlive the solver
and must not be modified after construction.

```cpp
using Solver = VectorSpace::block_wiedemann_solver<fieldType>;
Solver::options config; // block_size=8, trials=1, holdout=8, seed=17, threads=8
Solver solver(matrix, config);
auto result = solver.rank();
// result.rank, result.nullity, result.trial_ranks, result.probabilistic
```

The explicit operator constructor accepts dimensions and two callbacks implementing
multiplication by A and its transpose. Each callback receives a read-only input
span, an output span and a block width. Blocks are row-major arrays and callbacks
must overwrite the whole output. The native signed-byte adapter uses the field's
SmallSignedInt multiplication overloads.

## Rank algorithm

1. Work on the smaller side of a rectangular matrix, interchanging A and its
   transpose when needed. Independently sample nonzero diagonal scalings.
2. Apply B = D_left A D_right A^T implicitly. Ordinary A A^T is insufficient over
   finite fields: it can be zero even for a nonzero A.
3. Generate the full b-by-b sequence U^T B^(i+1) V from random dense blocks.
   Starting with B V removes the semisimple zero part of a generic preconditioned B.
4. Compute a shifted order basis of [S(z); I] using full block discrepancies.
   Extract a proper, row-reduced matrix recurrence; its determinantal degree is
   the sum of row degrees. Constant and leading coefficient matrices are checked
   for full rank. This replaces the old collection of scalar BM projections.
5. Validate the recurrence on held-out moments. If necessary, increase the
   training sequence up to a fixed dimension-based limit; failure is explicit.
6. Repeat with independent scalings and projections; return the largest observed
   rank and retain all trial ranks. Nullity always refers to A's original domain,
   even when the internal operator uses the transpose.

Rank results are **Monte Carlo estimates**, not deterministic certificates.
Agreement between trials and recurrence checks does not prove that projections
have observed every invariant factor. A recurrence failure throws; changing the
seed or block width is available through options. Small exact-elimination tests
are the correctness oracle, never an implicit fallback inside this solver.

The underlying block-sequence/preconditioning approach is described by
[Dumas, Elbaz-Vincent, Giorgi and Urbanska, PASCO 2007](https://www.lirmm.fr/~giorgi/pasco07.pdf).
This implementation uses an iterative order-basis algorithm and shifted sequence;
it does not implement their high-performance polynomial-matrix algorithms.

## Linear solves

`solve_MX_equals_y` keeps the existing unique_ptr-vector API. It includes the
scaled RHS in the starting block, obtains a right generator from transposed
moments, and reconstructs a candidate by Horner evaluation, with a single-vector
accumulator just as in the scalar solver. Every returned
solution is checked against the original rectangular A and y. An empty result
means no solution was found in the allotted trials, not proven inconsistency.
Zero RHS and empty dimensions are handled explicitly.

## Right nullspace

```cpp
auto kernel = solver.nullspace();
// kernel.basis: vector<vector<K>>, one domain-coordinate vector per column
// kernel.rank_estimate: the rank computation used to determine the target size
// kernel.complete: extracted the estimated nullity many independent vectors
```

This API works with the same signed-byte CSC, compressed sparse matrix, legacy
LIL and matrix-free constructors. For example, the callbacks can implement the
unweighted stack `[A; B^T]` and its transpose without allocating that stack.

Extraction uses the domain operator P = D_col A^T D_row A. Starting from random
V, a right block recurrence for P V reconstructs Z with P Z = P V; V-Z gives
candidate kernel vectors. Diagonal preconditioning is internal to Wiedemann;
the supplied operator is unchanged and output vectors use its original coordinates.
Every candidate is checked against the original A before admission. Incremental
echelon reduction removes dependencies and normalizes pivots, using the output
basis itself rather than retaining another dense copy.

For square/tall matrices, extraction reuses a right generator and starting block
from the rank trials. Wide matrices retain the smaller-side rank computation and
use a separate domain sequence for extraction. Reconstruction uses Horner
evaluation with only the number of still-needed output vectors, capped at the
block width. In addition to the
usual block workspace, output storage is O(domain_dimension * nullity); a large
nullspace can therefore be expensive even for a sparse input. Sampling is bounded
by ceil(estimated_nullity / block_size) + trials blocks. On exhaustion, the result
contains the verified independent vectors found so far and `complete` is false.
Rank recurrence failures retain the existing exception behavior. No field-retry
logic or homology interpretation is included.

`complete` is conditional on the Monte Carlo rank estimate, not a deterministic
completeness certificate. Each returned vector and its independence are checked
algebraically regardless of that estimate. A zero-row matrix returns the coordinate
basis; a zero-column matrix returns an empty basis.

## Cost and tests

Sparse data stay sparse; elimination fill-in in the original matrix is avoided.
Krylov blocks, moments and polynomial-generator storage grow roughly as dimension
times block width. The iterative order-basis phase has quadratic dependence on
sequence length. Recurrence state persists when extending a sequence, analogous
to scalar Berlekamp–Massey. Polynomial rows use flat storage, and each discrepancy
pivot inverse is computed once. Sparse products are parallel; recurrence and
projection work remain serial;
large-loop scaling still needs benchmarking and polynomial-stage optimization.

`tests/test_block_wiedemann.cpp` compares rank and nullity with exact sparse
elimination, checks nullspace basis size, independence and zero residuals,
and verifies solutions in 274 cases, covering empty, zero, singular,
rectangular, random low-rank matrices and an isotropic matrix with A A^T=0.
Block widths 1, 2 and 4 are tested; inconsistent systems and legacy LIL construction
are covered. Signed-byte CSC and an unweighted matrix-free stack are also tested.
Homology tests check ranks and nullspace bases of both contraction matrices
at every L6 stage in both parities against exact ranks.

The field concept in `Field.hpp` requires multiplication by `SmallSignedInt` in
both orders and in-place. All five finite-field implementations pass exhaustive
signed-byte multiplication tests, including ordinary-integer non-narrowing
regressions. Rational Q tests are conditional on Boost headers being available.

## Threads, transpose storage and profiling

Compile with `-fopenmp` (enabled by GraphHomology/Makefile). The sparse matrix
constructors use `options.threads`, default 8. Generic callback operators supply
their own threading; `GraphHomology::StackedContraction` accepts a worker count
as its third constructor argument, also default 8. Builds without OpenMP remain
serial. Products with fewer than 32,768 nonzero-times-block-lane operations run
serially to avoid thread-team overhead.

The GC stack stores transposes of down and up rather than assembling a separate
stacked matrix. Its transpose application accumulates both terms within the same
output-row iteration. No output scatter, atomic updates or per-thread dense
reduction arrays are used. The additional storage is one sparse transpose of each
differential: five bytes per nonzero plus column offsets.

`solver.timing` and `solver.print_timing()` expose accumulated Gram-product,
sequence/projection, recurrence and reconstruction timings, using the same
`timer_accum` utility as the scalar solver. Gram time overlaps the other phases;
do not add it to them. Detailed before/after logs and a repeated multiplication
benchmark are under `output/wiedemann_optimization`.

On the measured laptop, `OMP_WAIT_POLICY=PASSIVE` reduced overhead between sparse
products and the serial recurrence/projection work. It is an optional environment
setting, not forced by the solver; see the benchmark report for isolated timings.

Progress logging is optional: set `options.log = &std::cout` (or another stream
that outlives the solver). The graph matrix CLI enables it. Flushed lines report
rank trials, Krylov steps, recurrence construction and validation, Horner steps,
and verified independent vectors. Long phases report roughly every five seconds;
counters refer to the current phase, not a predicted overall completion time.
The default null stream keeps library calls quiet.

The default is one rank trial per run. To obtain independent runs, change the
seed: the graph CLI accepts `--nullspace OUTPUT_DIRECTORY [BLOCK_SIZE [SEED]]`
and `--rank [BLOCK_SIZE [SEED]]`. The seed controls both random diagonals and
projections and is recorded in the log and nullspace metadata. For example,
`--nullspace output/run_seed18 8 18` uses seed 18 with block size 8.

Recurrence discrepancy rows and recurrence validation samples also use
`options.threads` with OpenMP when there is enough work. Pivot selection and
dependent polynomial row elimination remain sequential. Validation uses bounded
batches so progress reporting stays on the controlling thread.

Krylov projections also use `options.threads`: the block moment entries are
independent dot products, with one writer per entry and no atomic updates or
shared reductions. Small projections stay serial to avoid OpenMP overhead.

For Z32783, Krylov projections use packed lane-major copies and delayed uint64
accumulation. The fixed projection is packed once per sequence, while a reusable
buffer receives the evolving block each step. Integer sums are reduced every
2^20 terms, well below uint64 overflow. The two packed buffers cost
`2 * gram_dimension * block_size * sizeof(Z32783)` bytes (about 104 MiB for L10).
Other fields retain the generic projection implementation.

Progress lines report both `elapsed_seconds` (steady-clock time, excluding
Linux system suspend) and `calendar_elapsed_seconds` (system-clock time,
including suspend). Both start at the same phase boundary. Calendar elapsed
time is affected by adjustments to the system clock; neither value is CPU time.
The five-second reporting cadence still uses the steady clock.


## Direct square operator

`block_wiedemann_solver<K>::from_square_operator(n, apply, options)` runs rank or
nullspace extraction directly on a supplied square map. It applies a random
nonzero left diagonal after each callback, preserving its right kernel, and
reuses the existing projected recurrence and Horner reconstruction. It does not
construct a transpose or an additional Gram product. The rank estimate assumes
that generic left diagonal scaling makes the zero eigenvalue semisimple, as for
the diagonally symmetrizable natural graph composition; this API is not a rank
guarantee for arbitrary square maps (e.g. strictly upper triangular maps).
Candidate vectors are checked against the original callback. The rectangular
linear-system solving method is unavailable in this mode.

The graph CLI's `--nullspace` now selects the transposed natural composition;
`--legacy-nullspace` retains the earlier stacked operator. Output conversion and
separate differential residual checks are documented in the graph README.

## Batch versus incremental recurrence

`options.incremental_recurrence` enables an opt-in streaming comparison; the
production default remains batch. Incremental mode feeds each new training term
into the existing block approximant-basis state while leaving the newest
`options.holdout` terms unprocessed. Candidate extraction is attempted every
`options.recurrence_check_interval` training terms (default 32), and also exactly
at the original batch boundary. This is interleaving on the controlling thread,
not concurrent CPU/GPU execution.

A structural candidate test is only a gate. Candidates first encounter the fresh
held-out terms, so bad predictions are rejected before rechecking all training
history. A successful candidate must still pass the original full validation,
including leading/constant block ranks. Failed checks continue from the retained
state. Nullspace reconstruction and original-operator residual checks are the
same in both modes; early stopping does not remove the probabilistic assumptions.

`sequence_stats` records sequence starts, generated moments, processed training
terms and candidate-validation attempts. Timing segments exclude recurrence
work from the Krylov multiplication/projection total even when interleaved.
The timing call count now counts those segments; use `sequence_stats.sequences`
to check whether reconstruction reused a sequence.

Reproduce the comparison on both L8,V12 parities and known-rank diagonal maps:

```bash
make -f GraphHomology/Makefile build/bench_block_recurrence HOMOLOGY_FLAGS='-std=c++23 -O3 -march=native -fopenmp -I. -IVectorSpace'
OMP_WAIT_POLICY=PASSIVE build/bench_block_recurrence output/triangle_comparison_L10_cut_splits_run1/splits_L8 5
```

Both modes use the same seed within a pair, block size 8, eight workers and eight
held-out terms. Execution order alternates between repetitions. Graph timing
includes nullspace extraction and its operator checks; basis/matrix construction,
coordinate conversion and canonical subspace comparison are outside the timed
region. Both converted differential residuals are verified, and canonicalized
output subspaces must match. Diagonal cases benchmark rank only and check the
known exact rank. No singular-solver code is changed.

## Packed storage and mmap archives

Projected moments, the working polynomial basis (recurrence bucket), and the
final generator use fixed contiguous coefficient slabs backed by `OwnedArray`.
Degree/length arrays are also contiguous. These buffers and recurrence scratch
are allocated before generating the first Krylov term. No polynomial row owns
another allocation. This guarantee covers recurrence storage, not every later
nullspace reconstruction allocation.

`options.sequence_capacity` sets the maximum number of projected blocks.
Zero reserves the initial batch plus one doubling, capped by the solver's
training limit, with holdout space added. At capacity exhaustion the solver
saves its state to `options.checkpoint_path` and throws, allowing the allocations
to be released. Restart with a larger capacity and `resume_checkpoint=true`.
The default algorithm still uses separate Krylov and recurrence stages.

`mmap_archive.hpp` provides a shared binary backend. Each component writes a
named, length-delimited section and has its own reader:

- Sparse CSC matrices: dimensions, offsets, graph indices and coefficients in
  separate typed arrays (`ContractionMatrixStorage::save/load`). Both int8
  contraction and int32 adjoint coefficients are supported.
- Vector blocks: initial vectors, projections, diagonals and current Krylov
  vectors use the same block serializer with distinct section names.
- Projected coefficients: used moment blocks, dimensions and field metadata.
- Recurrence state: progress, active polynomial lengths, shifted degrees and
  packed active coefficients. Unused reserved capacity is omitted.
- Final recurrence: row degrees and active coefficient blocks.

The archive header records format version, byte order and checksum algorithm.
Integer arrays record width, signedness and count. Each field block records
field name, characteristic, element width, dimensions and layout; each vector
lane also has its own field/width descriptor. Current solver blocks use one
field throughout: simultaneous mixed-prime computation is not implemented.
Solver checkpoints additionally record enumeration/coefficient/accumulator
widths, seed, solver dimensions, mode and progress counters.

Saving computes the exact file size, reserves disk space, maps a temporary file,
and bulk-copies the live contiguous arrays directly into it. Large payloads are
borrowed until `finish()`; keep their owners alive and unchanged until then.
After CRC32C, `msync` and `fsync`, the temporary file is renamed over the target.
Loading maps the file read-only, checks its checksum and metadata, then copies
payloads into the destination allocations. It is mmap-backed bulk loading,
not zero-copy mapped working storage. Native field representations require
compatible types and byte order; incompatible metadata is rejected.

Matrix archives are independently reusable; a solver checkpoint does not embed
an arbitrary operator callback. Recreate the same operator when resuming. The
loader compares regenerated initial vectors, projections and diagonals; this
is a consistency check, not a complete matrix identity proof. Only the current
full Krylov block is saved, together with all accumulated projected moments;
historical full Krylov vectors are not retained. Resume currently addresses
the first sequence of a run, using the same seed/options. Elapsed segment times
restart after loading.

For the graph executable, the checkpoint is `OUTPUT.recurrence`. Set
`GC_SEQUENCE_CAPACITY` for a custom capacity and `GC_RESUME_CHECKPOINT=1` when
restarting with the same output prefix. Sequence checkpoints are written on
capacity exhaustion, not periodically or on process termination. Reconstruction
has separate periodic checkpoints, described below.

Correctness checks (no benchmark):

```bash
make -f GraphHomology/Makefile build/packed_recurrence_test
build/packed_recurrence_test
```

These cover component round-trips, larger-capacity restore, checksum/type/field
rejection and matching uninterrupted versus resumed block-solver results.

## Reconstruction checkpoints

Nullspace reconstruction now has a separate archive at
`checkpoint_path + ".reconstruction"` (for the graph CLI:
`OUTPUT.recurrence.reconstruction`). It uses the same mmap backend and field/type
validation as the other components. `reconstruction_state.hpp` owns the
independently serializable Horner progress and accumulated-vector block.

The surrounding nullspace checkpoint contains the final recurrence, starting
vectors, initial operator image, reconstruction weights, diagonal scalings,
rank estimate, attempt number, already verified basis vectors and their pivots,
and the random-generator state. Temporary multiplication scratch and the
current coefficient combination are regenerated, not serialized.

A snapshot is written before the first Horner step, at complete-step boundaries
approximately every 60 seconds, and after the final step. Configure this with
`options.reconstruction_checkpoint_seconds`; zero saves every step. Snapshot
writing is synchronous. The period is checked between complete steps, so a
single long operator application can delay it. Each snapshot atomically replaces
the previous file. A hard stop may lose work since the last completed snapshot.

With `resume_checkpoint=true` (`GC_RESUME_CHECKPOINT=1` in the graph CLI), a
reconstruction archive takes priority over the earlier sequence checkpoint.
The rank/Krylov/recurrence stages are skipped. The operator must be recreated;
one saved starting-block image is recomputed as a consistency check. This is
not a full operator identity certificate. Horner resumes at the next saved
power, then the ordinary residual and independence checks run. Later attempts
continue with the saved random-generator state and accepted basis.

The final snapshot precedes verification of the current candidate block, so
reloading a completed reconstruction repeats those inexpensive checks, not
Horner evaluation. This integration is for `nullspace()`; the separate
single-right-hand-side solve API does not yet automatically checkpoint its
reconstruction. The existing running binary gains no checkpoint support until
it is rebuilt and restarted.


## Rank-to-nullspace handoff

With `options.checkpoint_path` set, `rank()` uses the domain-side operator and
reconstruction-compatible projection orientation. It saves the completed rank
result and best validated trial to `checkpoint_path + ".rank"`. The archive
contains rank/nullity, trial ranks, starting/projection blocks, initial operator
image, diagonal scalings and final recurrence, with solver/type/field metadata.
Ordinary rank calls without a checkpoint path retain the smaller-operator choice.
A checkpointed wide rectangular rank can therefore cost more than rank alone.

A fresh solver with the same operator/options and `resume_checkpoint=true`
loads this handoff in `nullspace()`, skipping rank estimation, Krylov generation
and recurrence finding for the first candidate block. It checks one saved
operator image before reconstructing. Additional candidate blocks may still
need new sequences. A reconstruction checkpoint takes precedence over a rank
handoff, which takes precedence over a sequence-capacity checkpoint.

The graph CLI separates preparation from extraction:

```bash
build/contraction_matrices_L9_V13 GRAPH_DIRECTORY even --prepare-nullspace OUTPUT 8 17
GC_RESUME_CHECKPOINT=1 build/contraction_matrices_L9_V13 GRAPH_DIRECTORY even --nullspace OUTPUT 8 17
```

Use the same `OUTPUT`, graph directory, parity, block size and seed. The handoff
is `OUTPUT.recurrence.rank`; the parent directory must exist. The preparation
command computes the rank of the actual natural-composition operator used for
representatives. The existing `--rank` command computes the two differential
ranks and does not provide this handoff. Rebuild before using the new option;
already running binaries are unchanged. Empty domain/codomain cases are handled
directly and do not need a saved recurrence.
