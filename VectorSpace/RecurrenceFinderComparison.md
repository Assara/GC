# Separate recurrence finders

## Algorithm reference

Claude-Pierre Jeannerod, Vincent Neiger, Gilles Villard,
*Fast computation of approximant bases in canonical form*.
Preprint: https://arxiv.org/abs/1801.04553
PDF: https://arxiv.org/pdf/1801.04553

The paper recalls PM-Basis, the divide-and-conquer approximant-basis algorithm.
This is the reference for the new experimental finder in
`divide_conquer_generator.hpp`. It is not an implementation of all the paper's
canonical-form algorithms, and does not claim their full complexity bounds.

For F = [S; I], approximation order K and row shift (0,...,0,1,...,1):

1. Compute P1 for the first floor(K/2) coefficients.
2. Form the residual from the next coefficients of P1 F, after removing the
   first floor(K/2) zero coefficients.
3. Compute P2 for that residual, using the shifted row degrees of P1.
4. Return P2 P1.

Small subproblems use iterative discrepancy elimination. Polynomial matrix
products use Karatsuba, with ordinary matrix coefficient multiplication at the
leaves. Operand order is preserved; matrix coefficients do not commute.
Residual multiplication currently computes a full product and extracts the
needed interval, rather than using a specialized middle-product kernel.

## Baseline and validation

`block_minimal_generator.hpp` remains the production reference finder. The
experimental finder imports its completed approximant basis into the existing
extraction and validation code. Full-history checks, held-out terms and the
constant/leading rank tests are unchanged. No validation shortcut is included.

Both finders use packed storage. Experimental scratch is allocated up front
and reused by recursion; memory is reported alongside runtime. It uses more
workspace than the reference implementation. Speedup is not assumed: compare
on the same sequence before choosing a production implementation.

## Capture and compare

Set `GC_RECURRENCE_SEQUENCE=PATH` on a newly built graph `--prepare-nullspace`
or `--nullspace` run to export its projected sequence before recurrence finding.
This is optional and adds no I/O when unset. A later training extension or
sequence replaces PATH with the latest input; use separate paths for runs you
want to retain. The archive records the training boundary, block size, total
length and field metadata, including all available held-out terms.

Build and run the standalone comparison after the active long run finishes:

```bash
make -f GraphHomology/Makefile build/bench_recurrence_finders HOMOLOGY_FLAGS='-std=c++23 -O3 -march=native -fopenmp -I. -IVectorSpace'
build/bench_recurrence_finders PATH 3 8 16 8 > recurrence_comparison.csv
```

Arguments are sequence file, repetitions, threads, PM-Basis leaf order and
Karatsuba leaf length. The same in-memory sequence is used for both finders.
Execution order alternates. Allocation/setup, recurrence construction and
unchanged validation are timed separately; workspace excludes the common input
and returned generator. Coefficient-by-coefficient equality and matching
acceptance are required. A recurrence rejection is reported in CSV; it is not a
successful extraction. Archive loading and Krylov generation are not timed.

The experimental finder is standalone and does not replace the production
solver or modify the singular solver. These two implementations deliberately
use matching deterministic pivot choices; mathematical validity alone does
not in general imply identical recurrence coefficients.

Correctness checks, without performance measurements:

```bash
make -f GraphHomology/Makefile build/recurrence_finders_test
build/recurrence_finders_test
```

The tests cover known-rank sequences, zero and singular cases, different block
sizes and recursion thresholds, odd orders, arbitrary sequences, corrupted
holdouts, and mmap input round-trips. No performance conclusion is recorded yet.

## Resumable divide-and-conquer construction

The finder now stores its work stack explicitly in a fixed owning array. Each
frame contains the approximation order, phase, leaf progress and offsets to its
input, output, degree arrays and partial results. There are no saved pointers.
All coefficient arenas and frame storage are allocated by the constructor;
`step()`, `process()` and loading into that object do not grow them.

`step()` performs one stable transition: frame setup, one leaf discrepancy
update, a residual multiplication, or a final basis multiplication. `process()`
continues from the current state rather than restarting. Its progress callback
runs at safe boundaries, where the caller can save and optionally throw to
pause. The number of processed leaf terms can reach the training length before
all basis multiplications finish; use `complete()` for completion.

`save(archive_output&)` and `load(archive_input&)` use the common mmap archive.
The `pm-basis-state` section includes format version, field name/characteristic
and element width, block/order/cutoff parameters, input length and identity hash,
progress, work-stack offsets, degree workspace, partial basis and live arena
prefix. Integer arrays record their widths. Stack topology and arena offsets
are validated when loading. The input sequence is supplied separately and must
match, including held-out terms; its FNV-1a identity hash is a consistency guard,
not a cryptographic certificate. The expanded input is regenerated from it.
Thread count may change on resume; approximation parameters must match.

Temporary Karatsuba/discrepancy scratch and unused arena capacity are not saved.
Snapshots are taken between full polynomial products, not during an inner
Karatsuba call. A hard stop loses work since the last completed snapshot; a
particularly long product can delay a periodic checkpoint. Validation still
runs after construction and is not checkpointed midway. A completed checkpoint
can be loaded and validated without repeating construction.

The standalone command periodically saves to `OUTPUT.state`, then writes the
validated packed recurrence to `OUTPUT`:

```bash
make -f GraphHomology/Makefile build/find_recurrence
build/find_recurrence SEQUENCE OUTPUT
build/find_recurrence SEQUENCE OUTPUT --resume
```

It writes a checkpoint at start, approximately every 60 seconds at a stable
boundary, and at construction completion. Save failures leave the previous
atomically published checkpoint available. The sequence and output parents
must exist. The comparison benchmark remains separate so checkpoint I/O is not
mixed into its recurrence-construction timing.

Resume tests destroy/recreate the finder at every stable transition for small
problems, verify equal step counts and identical output coefficients, test a
changed thread count, and reject changed source data. The standalone command
also passes a save/load round-trip on an exported GC sequence.
