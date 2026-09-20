# 32-bit bounds for the transposed natural composition

This concerns the proposed operator

    (S_down C_down + C_up S_up)^T

on graphs with V vertices and E edges. C is the integer contraction differential
and S is natural splitting, with consistent signs. It is NOT the currently
running diagonally preconditioned transpose-Gram implementation. No random
diagonals are included inside this proposed composition.

All graphs are simple and at least trivalent. Bounds count unsigned operation
paths before cancellations. Restricting to cut-free graphs or removing
orientation-zero graphs cannot increase the counts.

## Result for p=32783

Canonical residues have magnitude at most 32782. A signed 32-bit accumulator
can hold 2147483647, giving a sufficient absolute path-weight limit of 65508.

| Degree-3 target | V | E | Path-weight bound | Accumulator bound | Certified signed 32-bit |
|---|---:|---:|---:|---:|---|
| L10 | 14 | 23 | 3121 | 102312622 | yes |
| L11 | 15 | 25 | 7136 | 233932352 | yes |
| L12 | 16 | 27 | 15879 | 520545378 | yes |
| L13 | 17 | 29 | 34726 | 1138387732 | yes |
| L14 | 18 | 31 | 75069 | 2460911958 | no certificate from this bound |

For ALL feasible vertex numbers at a fixed loop order, the maximum of this
bound is 56826 at L11,V12,E22, still within 32 bits. At L12 it is 124315 at
V13,E24, outside 32 bits. Thus every simple at-least-trivalent graph sector
through L11 is covered, while individual sectors are covered further.

These are sufficient bounds. Failure of a bound does not prove overflow for
the actual matrix. The earlier conversational multiplication 3121*32782 was
mistyped; the correct value is 102312622, as in this table.

## Counting the composed columns

Set f(d)=2^(d-1)-d-1 and S(G)=sum_v f(d_v). Splitting first and contracting
second has at most (E+1) S(G) paths.

For contraction of an edge with endpoint valences a,b, a valid simple target
has split count

    S(G) - f(a) - f(b) + f(a+b-2).

Summing this expression over all source edges upper-bounds contraction followed
by splitting. Invalid contractions count as zero; including their positive
hypothetical contributions only enlarges the bound.

Define

    g(a,b) = f(a+b-2)-f(a)-f(b).

Then g(3,3)=3 and g(d,3)=2^(d-1)-1. Write every nontrivalent valence as
3+t_i, with t_i>=1, and put u_i=2^t_i-1. The extra contribution from an edge
joining two nontrivalent vertices, beyond accounting for its endpoint degrees
separately against trivalent neighbors, is exactly

    g(3+t_i,3+t_j)-g(3+t_i,3)-g(3+t_j,3)+3 = 8 u_i u_j.

By simplicity, each such vertex pair supports at most one edge. Therefore

    contract_then_split <= E(S(G)+3)
                           + 4 sum_i (t_i+3)u_i
                           + 8 sum_{i<j} u_i u_j.

Adding (E+1)S(G) bounds the sum of the two natural compositions. Here

    S(G) = 4 sum_i u_i - sum_i t_i,
    sum_i t_i = 2E-3V,
    1 <= t_i <= V-4,
    number of t_i <= V.

Maximizing this explicit expression over partitions of 2E-3V, with those
constraints, gives a bound depending only on V,E. Some partitions may not be
graphical, and some pairs cannot simultaneously support edges; retaining them
makes the result conservative. For L10,V14,E23 the maximum is attained in this
relaxation by one excess of four, i.e. valences (7,3,...,3), giving 3121.

## Why transposing certifies output accumulators

For each untransposed composition, the operation-path column bound dominates
the corresponding column sum of the product of the entrywise absolute factor
matrices. Transposing turns that into an absolute path ROW bound. Multiplying
by the maximum input magnitude therefore bounds final accumulators and partial
sums in the SECOND stage, regardless of summation order. Summing the two
nonnegative bounds also covers accumulation of both branches before reduction.

First-stage intermediates can separately be bounded by the column bounds of
natural splitting on the lower space and contraction on the upper space:
Smax(V-1,E-1) and E+1. These also fit signed 32 bits in the certified cases
above and are checked in the certificate script. Any intermediate that does
not connect to a valid second-stage output can instead be omitted.

This reasoning requires signed accumulation, promotion before multiplication,
and normalization of negative remainders at the final reduction. It does not
permit inserting unreduced random diagonal multipliers without updating the
bounds. It also does not prove that this proposed natural operator has the
same right kernel as the currently implemented stacked operator.

## Reproduction

```bash
python3 tools/check_natural_composition_bounds.py --prime 32783 --max-loop 14
```

The script enumerates degree-excess partitions, not graphs, and uses exact
integers. Full results are saved in
`output/accumulator_bounds/natural_composition_prime32783.json`.

## Implementation

`NaturalAdjoint.hpp` constructs the adjoint from the contraction matrix and
integer automorphism sizes. For C: source -> target, let D_source and D_target
be diagonal matrices of automorphism group orders. Then

    S = D_source^{-1} C^T D_target
    S[G,H] = C[H,G] * |Aut(H)| / |Aut(G)|.

This is the adjoint for the diagonal pairing <G,G> = |Aut(G)|. Counting an
edge contraction together with an isomorphism of its target gives the same
objects as counting a vertex split together with an isomorphism of its source:

    signed_contraction_count * |Aut(H)|
      = signed_split_count * |Aut(G)|.

The existing convention appends the new edge and vertex at the end when
splitting. Contracting that new edge has sign +1, so there is no additional
degree-dependent global sign. Tests compare every adjoint column against the
existing split differential on L6, both parities, with cut-vertex and
orientation-zero graphs projected out.

Sizes come from `standardize4_with_aut_count`, once per basis graph. The older
`automorphism_group_size` routine is not used: on L6 even V8 basis index 6 it
reports 2 instead of 1. The newer counts are independently checked by enumerating
all vertex permutations for the L6 V8 bases in both parities.

Construction cancels the automorphism ratio using gcd, checks divisibility,
then performs a checked integer multiplication. There is no division in the
finite field and no production split enumeration. The CSC sparsity pattern is
exactly that of the contraction transpose. Contraction coefficients stay int8;
natural adjoint coefficients use int32 because the contraction edge-count bound
does not apply to them. All final arrays use `OwnedArray`.

`NaturalComposition<K, Accumulator>` evaluates the transposed sum above, using
row-owned OpenMP gathers with reusable integer workspaces. Default accumulator:
`GraphAccumulator = std::int32_t` in `types.hpp`. Change that alias to
`std::int64_t`, or explicitly instantiate `NaturalComposition<K,std::int64_t>`.
Only the final output conversion reduces modulo p. No random diagonals are
inserted inside this operator.

At construction, the evaluator computes actual absolute path-row bounds from
the four stored factors, including both branches and the first-stage column
sums. For each middle column c, the second-stage bound is

    sum_h |C_down[h,c]| * sum_g |S_down[g,h]|
      + sum_u |S_up[u,c]| * sum_g |C_up[g,u]|.

Multiplication by p-1 must fit the selected signed accumulator. Bounds are
computed with saturating arithmetic, so even a rejected certificate cannot
overflow. An insufficient int32 certificate throws with instructions to use
int64. There are no overflow checks in the multiplication loop. Each evaluator
borrows the matrices and owns mutable scratch space; it is not reentrant.

The CLI option `--natural-adjoint [BLOCK_SIZE [SEED]]` builds these adjoints,
reports storage and the actual accumulator certificate, and evaluates one
sample block. `--rank` retains the original differential rank calculation. `--nullspace` now
uses this transposed composition directly, converts by inverse automorphism
orders and verifies both contraction and natural splitting residuals. The old
stack path is available as `--legacy-nullspace`.

Example (L6,V8):

```bash
make -f GraphHomology/Makefile BINARY=build/contraction_matrices_L6_V8_adjoint
build/contraction_matrices_L6_V8_adjoint output/triangle_comparison_L10_cut_splits_run1/splits_L6 even --natural-adjoint
```

### Integrality verification

Every nonzero entry is checked during construction, in all builds (not just
assert-enabled tests). With a=|Aut(G)|, b=|Aut(H)| and g=gcd(a,b), the necessary
and sufficient condition is

    (a/g) divides C[H,G].

The reduced numerator b/g is coprime to a/g. Therefore the remainder check on
|C[H,G]| before division exactly tests integrality of C[H,G]*b/a. Zero entries
are automatically integral. A failed check throws with source/target indices,
the contraction coefficient, and both automorphism sizes; no fractional value
is rounded or truncated.

The signed counting identity above proves this divisibility for the natural
adjoint on the retained oriented graph bases, assuming correct contraction
coefficients and automorphism orders. Empirical tests additionally check these
implementation assumptions.

`tools/check_natural_adjoint.cpp` independently checks every nonzero via

    S[G,H] * |Aut(G)| == C[H,G] * |Aut(H)|

using signed 128-bit cross-products, with no division or modular arithmetic.
It checks both adjacent matrices and both parities. For example:

```bash
make -f GraphHomology/Makefile HOMOLOGY_LOOP=9 HOMOLOGY_VERTICES=13 build/check_natural_adjoint_L9_V13
build/check_natural_adjoint_L9_V13 output/triangle_comparison_L10_cut_splits_run1/splits_L9
```

### Solver integration and output coordinates

The block solver's `from_square_operator` factory evaluates this composition
once per Krylov step, without forming another Gram product. Its random nonzero
left diagonal is applied in the field after the composition has reduced. This
preserves the right kernel; the rank estimate retains the generic semisimple-zero
assumption. Candidate residuals are tested against the supplied unscaled operator.

With D_V=diag(|Aut(G)|), the transposed splitting differential is
`S_down^T = D_(V-1) C_down D_V^-1`. Thus `x=D_V^-1 z` changes from the transposed
complex back to the original contraction complex, including boundaries. All D
entries must be nonzero in the field. After conversion, the CLI checks both
`C_down x=0` and `S_up x=0`, then normalizes the first nonzero coefficient.
A kernel of the sum alone need not satisfy the separate residuals over a finite
field; failing vectors cause an error, not an exported class.
