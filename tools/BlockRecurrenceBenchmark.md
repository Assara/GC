# Batch and incremental block recurrence comparison

Both implementations use the same field F32783, seed, block size 8, eight
workers, and eight held-out terms. Five seed pairs (17–21) were run sequentially,
with batch/incremental execution order alternating. The build uses
`-O3 -march=native -fopenmp` and `OMP_WAIT_POLICY=PASSIVE`. Incremental updates
happen after every new training term; stopping checks happen every 32 terms and
at the original batch boundary. This is CPU interleaving, not concurrent GPU/CPU
execution. Other machine activity was not controlled; small timing differences
should not be interpreted as a reliable speedup.

| Case | Rank / dimension | Krylov terms batch → incremental | Median batch | Median incremental |
|---|---:|---:|---:|---:|
| Diagonal, low rank | 64 / 2048 | 522 → 40 | 0.1151 s | 0.0065 s |
| Diagonal, nearly full rank | 2047 / 2048 | 522 → 520 | 0.3510 s | 0.3802 s |
| L8,V12 even | 2893 / 2894 | 734 → 734 | 1.5009 s | 1.4666 s |
| L8,V12 odd | 2639 / 2639 | 670 → 670 | 1.1443 s | 1.1401 s |

Diagonal timings measure rank computation only; graph timings measure complete
nullspace extraction, excluding graph/matrix construction and subsequent graph
coordinate conversion. The even graph case has nullity one, and the odd case
has nullity zero. Converted representatives pass both differential residuals;
canonicalized output subspaces agree for every seed pair. All diagonal ranks
agree with their known exact values. The block-solver regression suite passed
299 cases, including streaming modes, recurrence reuse and corrupted holdout
rejection. SHA256 checks confirmed the singular solver sources were unchanged.

The term count is the robust result: early stopping saves substantial work for
a genuinely low-rank map, but none on these graph examples. Timings are comparable
on the graph cases. Batch remains the default; incremental mode is available
through `options.incremental_recurrence=true`.

An initial incremental implementation rescanned training history before trying
the held-out terms on each candidate; that made the graph cases substantially
slower. The retained implementation tests fresh terms first and only performs
full validation on candidates that survive. It does not weaken acceptance:
accepted candidates still pass the original full-history and block-rank checks.

Raw results and configuration:

- `output/block_recurrence_holdout_first_9zgvnmo1/results.csv`
- `output/block_recurrence_holdout_first_9zgvnmo1/summary.json`
- `output/block_recurrence_holdout_first_9zgvnmo1/config.json`
- Initial full-history-first comparison: `output/block_recurrence_comparison_8seaca7g/results.csv`

Reproduction:

```bash
make -f GraphHomology/Makefile build/bench_block_recurrence HOMOLOGY_FLAGS='-std=c++23 -O3 -march=native -fopenmp -I. -IVectorSpace'
OMP_WAIT_POLICY=PASSIVE build/bench_block_recurrence output/triangle_comparison_L10_cut_splits_run1/splits_L8 5
```
