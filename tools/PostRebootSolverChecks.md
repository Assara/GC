# Post-reboot solver checks and benchmarks

Revision tested: `094df12` (no solver source changes during these checks).
CPU: Intel Core i7-1165G7, 4 cores / 8 hardware threads, 16 GiB RAM.
Build: C++23, GCC, `-O3 -march=native -fopenmp`, assertions enabled in tests.
Runs use `OMP_WAIT_POLICY=PASSIVE`. Benchmarks execute sequentially; host activity
and CPU frequency are not controlled. Medians use three alternating-order pairs.

## Correctness

- Signed-byte field multiplication passes for all five supported prime fields.
  Rational-field checks were skipped because Boost headers are unavailable.
- 299 block-Wiedemann rank/nullspace/verified-solution cases pass.
- 98 recurrence comparison and checkpoint-resume cases pass, including
  recreating the divide-and-conquer finder at every stable transition.
- Packed storage, mmap serialization, metadata/checksum rejection, capacity
  restart, rank handoff and interrupted reconstruction tests pass.
- L6 graph differential, integer d-squared, direct-splitting adjoint, 32/64-bit
  composition, overflow and representative checks pass.
- L8,V12 rank handoff works on both parities. Even parity gives one verified
  representative; odd parity has nullity zero. Resumed extraction performs no
  new Krylov generation or recurrence finding.

## Interpretation

Short graph sequences favor the reference finder. Longer synthetic sequences
show a modest divide-and-conquer advantage at substantially higher memory cost.
These are not L10 runtime estimates. Synthetic inputs are dense matrix sequences
obtained by left/right invertible mixing of eight scalar sums of distinct
nonzero geometric sequences over F32783. Their coefficient matrices are dense,
but their simultaneous diagonal structure is simpler than a general graph case.
The generator source is archived alongside the measurements for reproduction.

Validation is unchanged and measured separately. Every comparison requires
identical coefficients and matching acceptance. Workspace excludes the common
projected sequence and returned generator. No speedup from reduced validation
is included.

The L10-shaped projection benchmark measures 1,703,974 rows and width 8. Packed
uint64 remains the fastest tested projection kernel. Its reported kernel time
excludes packing: the benchmark reports packing both input blocks separately,
whereas production keeps the fixed projection block packed and repacks only
its changing block. These numbers are not complete Krylov-step timings.

## Before a long restart

Periodic Krylov checkpoints are still missing from the production path; it
currently saves its sequence state on capacity exhaustion. The experimental
PM-Basis finder has periodic checkpoints through its standalone driver but is
not selected automatically by the graph solver. Those integration points
remain separate from the passing checkpoint component tests.

The big computation was not restarted.

## Recorded medians (eight threads)

| Input | Training terms | Reference construction | PM-Basis construction | Reference / PM workspace |
|---|---:|---:|---:|---:|
| L8,V12 even | 726 | 0.0712 s | 0.1589 s | 0.72 / 7.66 MiB |
| L8,V12 odd | 662 | 0.0612 s | 0.1555 s | 0.66 / 7.42 MiB |
| Synthetic degree 512 | 1026 | 0.1089 s | 0.3373 s | 1.01 / 13.75 MiB |
| Synthetic degree 2048 | 4098 | 4.4579 s | 3.5144 s | 4.01 / 55.00 MiB |
| Synthetic degree 4096 | 8194 | 19.5998 s | 14.7534 s | 8.01 / 110.00 MiB |

L10-shaped packed uint64 projection median: 12.69 ms with eight threads,
30.74 ms with one thread, excluding packing.

Raw logs and CSVs: `output/solver_checks_20260920_202627`.

Decision: set divide-and-conquer aside; keep the reference finder in production.
No further PM-Basis benchmarking or integration is planned for this restart.
