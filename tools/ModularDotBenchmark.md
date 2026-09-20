# Exact modular dot-product experiment

This standalone benchmark specializes kernels to F32783. Running it does not
modify a solver process. The packed integer kernel is now shared with the
production Krylov projections in `VectorSpace/packed_projection.hpp`; the
floating-point alternatives remain experimental.

```bash
make -f tools/Makefile.modular_dot test
make -f tools/Makefile.modular_dot bench
```

The default flags are `-O3 -march=native -fopenmp`, without `-ffast-math`.
To compare the solver's existing compilation flags:

```bash
make -f tools/Makefile.modular_dot DOT_BINARY=build/bench_modular_dot_O1 \
  DOT_FLAGS='-std=c++23 -O1 -fopenmp -I.'
OMP_WAIT_POLICY=PASSIVE build/bench_modular_dot_O1 --case 1703974 8 8
```

`--case N WIDTH THREADS` benchmarks one shape. The default run measures L9 and
L10 projection shapes with one/eight threads, then a contiguous dot with
1,703,976 terms (approximately the L10 recurrence degree times block width).
That last case measures the arithmetic kernel, not the full recurrence algorithm.
`--test-only` performs correctness checks; `--packed-only` times packed layouts.

## Kernels

- `field`: original field multiplication and addition, reduction at every term.
- `uint64`: integer accumulation, modular reduction after each chunk.
- `double_cast_uint64_mod`: double accumulation, convert to uint64 for reduction.
- `double_reciprocal`: double accumulation, reciprocal quotient estimate,
  subtract the multiple of the prime, correct the remainder.
- `preconverted_double`: same reciprocal kernel with double inputs. Still strided.
- `packed_uint64`: integer accumulation on lane-major packed field inputs.
- `packed_double`: double accumulation on lane-major packed double inputs.

Canonical inputs lie in [0,32782]. All accumulation kernels reduce at most every
2^20 terms. The maximum partial sum is below 2^51, so every integer product and
partial sum is exactly representable in double. At this bound the floating
quotient estimate differs by at most one; corrections enforce [0,32783).
The generic Field implementation remains untouched. These specialized kernels
must not be reused with a different prime without deriving new bounds.

The packed layout makes each projection dot product contiguous. It is additional
storage, not a replacement for the row-major sparse-product blocks. Both-input
packing/conversion time is reported separately and is excluded from kernel times;
allocation time is also excluded. A real integration must account for conversions
per step, buffer allocation, and the ability to cache the fixed projection block.

## Checks and measurement

Checks cover random values, all-maximum residues, empty input, unequal strides,
chunk boundaries and multiple chunks, and more than 300,000 reciprocal-reduction
checks around random multiples of the prime and near the supported upper bound.
Each timed projection is compared exactly with the field-arithmetic reference.

Each shape uses one warm-up round and five measured rounds, rotating method order.
Output reports median, minimum, and maximum wall time. `OMP_WAIT_POLICY=PASSIVE`
is used consistently. Floating kernels use OpenMP SIMD reductions; exact bounds
make reassociation safe without relying on approximate arithmetic.

Initial logs are in `output/modular_dot_benchmark/`. They were collected while the
user's L10 solve was running, so CPU/cache/bandwidth contention and clock changes
limit conclusions. They are kernel timings, not end-to-end solver speedups.
