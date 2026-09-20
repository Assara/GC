# Exact accumulator bounds for GC contraction matrices

The proof below bounds integer accumulation before modular reduction. It uses
**the signed-unit edge contractions and inverse vertex splits**, not the bit width
used to store matrix coefficients. It applies to either parity convention.

The conclusions concern the current simple-graph complex, with minimum valence
three and parallel-edge contractions discarded. Removing cut-vertex graphs or
orientation-zero graphs cannot increase the bounds. A multigraph complex would
require a different inverse-split argument.

## Results

Let **E be the edge count in the middle basis**. The two contraction matrices
in the window involve graphs with E-1, E, and E+1 edges. Thus, for the running
L10,V14 case, E=23 and the largest contraction source has 24 edges.

Let `A = [down; up^T]`. The running solver applies

    x -> D_col A^T D_row A x,

where both diagonal matrices have canonical nonzero residues modulo p.
All limits below refer to **signed** 32-bit or 64-bit accumulators and canonical
input residues in [0,p-1]. They certify every middle edge count from 6 through
the stated limit. They are sufficient bounds, not claims of optimality: the
first edge count outside a limit is where this proof becomes inconclusive,
not where an actual matrix necessarily overflows.

| Reduction schedule | Every p < 2^16, signed 32 | Every p < 2^16, signed 64 | Every p < 2^32, signed 32 | Every p < 2^32, signed 64 |
|---|---:|---:|---:|---:|
| Reduce after each sparse product, before its diagonal | E <= 23 | E <= 83 | — | E <= 53 |
| Fuse each sparse product with its following diagonal, then reduce | — | E <= 53 | — | — |
| Carry both sparse products and the middle diagonal; reduce before the outer diagonal | — | **E <= 24** | — | — |
| Carry both products and both diagonals before any reduction | — | E <= 10 | — | — |
| For comparison only: carry both products with no diagonals | E <= 10 | E <= 39 | — | E <= 24 |

A dash means no uniform certificate for nontrivial simple-GC edge counts from
this bound. In particular, canonical residues for primes near 2^32 do not even
fit in a signed 32-bit accumulator. For such primes, diagonal multiplication
can exceed signed 64 bits even after both operands have been reduced. It fits
in **unsigned** 64 bits: (p-1)^2 < 2^64. Use that or a wider scalar multiplication
for the separately reduced diagonals in the first row of the table. This does
not make unsigned accumulation appropriate for signed sparse sums.

The comparison row without diagonals is **not** a proposal to remove the
solver's random preconditioning.

For our particular p=32783, the same certificates give:

| Reduction schedule | Signed 32 | Signed 64 |
|---|---:|---:|
| Each sparse product, before its diagonal | E <= 25 | E <= 85 |
| Each sparse product fused with its diagonal | — | E <= 57 |
| Both products and the middle diagonal | — | **E <= 26** |
| Both products and both diagonals | — | E <= 13 |
| Both products without diagonals | E <= 11 | E <= 40 |

## 1. Column bound from edge contraction

Write D_m for the integer contraction matrix from graphs with m edges to graphs
with m-1 edges. The implementation uses one oriented graph per isomorphism
class, with no division or multiplication by automorphism-group orders.

Each of the m source edges contributes either zero or one signed target graph.
Canonicalization changes only the sign; merging identical targets and cancelling
opposite signs can only decrease the sum of absolute coefficients. Therefore,
for every column G,

    sum_H |(D_m)_{H,G}| <= m.                         (1)

This also proves each individual integer coefficient has absolute value at
most m. It makes no assumption that an aggregated coefficient is still +/-1.

The row sums of D_m are different: many distinct source graphs may contract to
one target. Equation (1) alone does not bound a row of D_m by m.

## 2. Inverse-split count

Fix a target graph H and a vertex of degree d. Splitting that vertex adds an
edge between two new vertices. Each new vertex must receive at least two of
the d old incident edges, to remain at least trivalent. The two parts are
unordered, so the number of possible partitions is

    f(d) = (2^d - 2 - 2d) / 2 = 2^(d-1) - d - 1.     (2)

In particular f(3)=0. Every nonzero parent G in a row of D_m arises from one
such split of H: choose an isomorphism of the contracted graph with H, and
undo the contracted edge. There may be many splits yielding the same parent;
there cannot be more distinct parents than splits.

For S(H) = sum_v f(deg(v)), there are at most S(H) distinct parent columns.
Each has coefficient at most m by (1). Hence

    sum_G |(D_m)_{H,G}| <= m S(H).                    (3)

This deliberately conservative factor m accounts for multiple source edges
contracting to the same unlabeled target, including automorphism-related
multiplicities. No assumption about asymmetry is needed.

## 3. Bound S(H) using only its edge count

Suppose H has e edges and n vertices. Simplicity and minimum degree three give

    4 <= n <= floor(2e/3),
    3n <= 2e <= n(n-1),
    3 <= d_v <= n-1,       sum_v d_v = 2e.             (4)

Relax the degree list by dropping graphicality and connectedness. This enlarges
the set over which we maximize S(H), so still gives an upper bound.

The increment f(d+1)-f(d)=2^(d-1)-1 is increasing. If two degrees lie strictly
between 3 and n-1, move one degree unit from the smaller to the larger. The
sum of f cannot decrease. Repeating this operation shows a maximizing list has
all but at most one degree at one of the endpoints.

For n>4 define

    t = 2e - 3n,
    a = floor(t/(n-4)),
    r = t - a(n-4),       0 <= r < n-4.

The maximum of the relaxed problem is

    F(e,n) = a f(n-1) + f(3+r),                       (5)

where the last term is zero for r=0 and is omitted if a=n. These conventions
agree because f(3)=0. For n=4 the only possible graph is K4, and F(6,4)=0.

Let S(e) be the maximum of F(e,n) over the feasible n in (4), or zero if there
are none. Define the monotone envelope

    B(e) = max_{0 <= j <= e} S(j).                    (6)

Then every eligible graph with at most e edges satisfies S(H) <= B(e).
Equations (2), (4), (5), and (6) are an explicit finite computation in integers;
no graph enumeration or experimental estimate is required.

## 4. Bounds for the stacked operator

In A=[down;up^T], down=D_E and up=D_(E+1). Combining (1) and (3) gives

    ||A||_infinity <= R(E) = max(E+1, E B(E-1)),       (7)
    ||A^T||_infinity <= T(E) = E + (E+1) B(E).        (8)

Here ||M||_infinity is the largest absolute row sum. For A, rows come either
from down (bounded by E B(E-1)) or up^T (bounded by E+1). For A^T, each output
row combines a down column and an up row, so their bounds must be ADDED.

Set H(E)=max(R(E),T(E)) and J(E)=R(E)T(E). All four functions are nondecreasing
because B is nondecreasing. This is why checking successive edge counts and
stopping at the first failed bound certifies every smaller edge count.

## 5. General prime p and accumulator width w

Set q=p-1 and M_w=2^(w-1)-1. Signed input contributions can be negative, so we
bound their absolute sums. If that bound is at most M_w, every partial sum in
any order is also within [-M_w,M_w]. The bound also covers each individual
matrix-coefficient multiplication. We do not rely on cancellation.

| Schedule | Sufficient bound on the absolute accumulator value |
|---|---:|
| One sparse product, then reduction | q H(E) |
| One sparse product followed by its diagonal, then reduction | q^2 H(E) |
| Both products, no diagonals | q J(E) |
| Both products with the middle diagonal | q^2 J(E) |
| Both products with both diagonals | q^3 J(E) |

Proof: canonical inputs have absolute value at most q. The first product is
bounded by q R. Its diagonal multiplies that by at most q, giving q^2 R.
The second sparse product multiplies the absolute bound by at most T, giving
q^2 R T. The final diagonal gives q^3 R T. Reducing between the products resets
the input bound to q, which gives the first two rows. Omitting diagonals gives
the comparison row. The same product J bounds the reversed AA^T order.

For a row of this table with bound q^k C(E), the general certificate is

    (p-1)^k C(E) <= 2^(w-1)-1,                       (9)

or equivalently

    p <= 1 + floor( (floor(M_w/C(E)))^(1/k) ).        (10)

The root in (10) is an integer root, not a floating-point approximation.
For any prime below the right-hand side, the certificate holds. Primality is
not needed for the overflow proof itself; it is an assumption of the solver.

Separately reduced diagonal multiplications require a type supporting q^2,
even when only q H(E) is used to certify the sparse accumulators.

For a uniform prime-bit bound p<2^b, use q <= 2^b-2 in (9). This is deliberately
conservative and avoids relying on which number just below 2^b is prime. The
16-bit and 32-bit tables above use b=16 and b=32 respectively.

## 6. Numerical certificate for E=23, p=32783

The integer formulas give

    B(22) = 1012,       B(23) = 1022,
    R(23) = 23276,      T(23) = 24551,
    q = 32782.

Consequently:

| Schedule | Absolute accumulator bound |
|---|---:|
| Each sparse product | 804,830,882 |
| Each sparse product fused with its diagonal | 26,383,965,973,724 |
| Both products, no diagonals | 18,733,243,609,432 |
| Both products with the middle diagonal | 614,113,192,004,399,824 |
| Both products with both diagonals | 20,131,858,660,288,235,030,368 |

Compare these with

    INT32_MAX = 2,147,483,647,
    INT64_MAX = 9,223,372,036,854,775,807.

Thus signed 32-bit accumulation is certified for each sparse product separately.
Signed 64-bit accumulation is certified through BOTH sparse products AND the
middle diagonal, provided we reduce BEFORE multiplying by the final diagonal.
The latter unreduced multiplication is not certified in signed 64 bits.
After reducing to [0,q], the final scalar multiplication is at most
q^2=1,074,659,524, which does fit in signed 32 bits for this particular prime.

At E=23 the upper bound on p from (10) for carrying the middle diagonal through
both products with signed 64-bit accumulators is p<=127045. Thus this schedule
is safe for every prime below 2^16, not just our current prime.

## 7. Reproduce and apply

Run from the repository root:

```bash
python3 tools/check_contraction_accumulator_bounds.py --prime-bits 16 --edges 23
python3 tools/check_contraction_accumulator_bounds.py --prime-bits 32 --edges 23
python3 tools/check_contraction_accumulator_bounds.py --prime 32783 --edges 23
```

The script uses arbitrary-precision integer arithmetic and integer roots. It
also independently checks (5) by dynamic programming over all bounded degree
lists for n=4 through n=10. It prints the bounds, signed limits, maximal prime
bounds at the requested E, and the first-threshold edge certificates. In
prime-bit mode `modulus_bound` is the upper bound 2^b-1, not a claim that this
number is prime.

Implementation requirements:

- Promote operands BEFORE multiplication; assigning an overflowed narrow product
  into a wide accumulator does not fix it.
- An unreduced intermediate must use a signed integer buffer, not a field element
  whose invariant requires values in [0,p-1].
- Normalize negative remainders back into [0,p-1] at each stated reduction point.
- Retain random diagonals. Their factors of q are included in this proof.
- Preserve the integer contraction coefficients exactly. This document does not
  use or certify any particular coefficient storage format.

This document and certificate do not modify the multiplication kernels or any
running computation. Actual absolute row sums of a constructed matrix can give
much tighter certificates than this uniform graph-combinatorial bound.

## Concrete reason split count alone is insufficient

The even L6 contraction matrix from V10,E15 to V9,E14 contains a row whose
absolute coefficient sum is 15, although the target has degrees
(4,3,3,3,3,3,3,3,3) and therefore only three candidate vertex splits.
Only one surviving parent basis graph contributes to that row, with coefficient
+15: all fifteen source-edge contractions produce the same oriented target.

In the saved sorted bases this is source column 5, target row 12. The source
edge list is:

    01 02 03 14 15 26 27 38 39 47 49 56 58 69 78

Reducing the number of distinct parents therefore does NOT compensate for
contraction multiplicities in the unweighted isomorphism-class basis. Splits
of a fixed target and contractions of a fixed parent count different objects;
a one-to-one counting argument needs automorphism weights, which our matrix
does not use.

Reproduce with `tools/check_gc_contraction_multiplicity.cpp` and the generated
L6 inputs. Observed output is saved in
`output/accumulator_bounds/multiplicity_L6.log`. This counterexample is not used
to establish the general bounds; it explains why replacing (3) by S(H) alone
would be invalid.

## Sharper fixed-(V,E) bound using split endpoint valences

The uniform E-times-split-count bound above can be improved without knowing the
actual graph or its automorphism group. This refinement applies to contraction
from source (V,E) to target (V-1,E-1).

Let h_d be the number of target vertices of valence d. A split at a degree-d
vertex partitions its incident edges into k and d-k, yielding endpoint valences

    a = k+1,   b = d-k+1,   a+b-2 = d.

Any edge of that parent contracting back to the target must have endpoint
valences {a,b}. To prove this, compare the parent and target valence histograms:
the contraction removes valences a,b and adds d, with d strictly larger than
either endpoint valence because both are at least three. The largest changed
valence therefore uniquely identifies d, and the remaining changes identify
the multiset {a,b}.

Since a,b<d, for a!=b the parent contains N_a=h_a+1 and N_b=h_b+1 vertices of
those valences. The number of eligible edges is at most

    M(h,a,b) = min(E, N_a*N_b, a*N_a, b*N_b).

For a=b, the parent has N_a=h_a+2 such vertices, giving

    M(h,a,a) = min(E, N_a*(N_a-1)/2, floor(a*N_a/2)).

These follow respectively from simplicity and the available incident half-edges.
They replace the blanket multiplicity bound E by a bound specific to that split.
For example, if a and b are distinct and neither occurs in the target, each
occurs exactly once in the parent: only ONE edge could have that endpoint type.

For 2<=k<=floor(d/2), the number of unordered partitions of type k,d-k is

    w(d,k) = binomial(d,k)       if 2k<d,
             binomial(d,k)/2     if 2k=d.

Consequently a target row has absolute sum at most

    F(h) = sum_d h_d * sum_k w(d,k) * M(h,k+1,d-k+1).

As before, multiple splits may yield the same parent. Counting every split
separately only overestimates the row sum. This argument retains multiplicities;
it does not assume they cancel or divide them by automorphism orders.

A rigorous bound depending only on V,E is obtained by maximizing F(h) over

    sum_d h_d = V-1,
    sum_d d*h_d = 2(E-1),
    3 <= d <= V-2,  h_d nonnegative integers.

Enumerating partitions of the excess 2E-3V+1 enumerates exactly these degree
multisets. Including non-graphical multisets makes the maximization conservative.
This is still not asserted to be optimal.

For source V14,E23 there are seven such multisets. Their maximum F(h) is 304,
at h_5=1, h_4=3, h_3=9:

- The degree-five vertex has 10 splits of type (3,4), each bounded by 16
  matching-valence contraction edges: contribution 160.
- Each of three degree-four vertices has 3 splits of type (3,3), each bounded
  by 16 matching-valence contraction edges: contribution 144.

This is NOT the degree distribution maximizing the unweighted split count.
It is the maximum of the refined multiplicity-weighted upper bound.

For source V15,E24 the corresponding maximum is 234. Hence the running stack
has the certified bounds

    ||A||_infinity <= max(304,24) = 304,
    ||A^T||_infinity <= 23+234 = 257.

At p=32783, separate sparse-product accumulators are therefore bounded by

    304 * 32782 = 9,965,728.

This improves the earlier fixed-(V,E) bound of 89,722,334 and the edge-only
bound of 804,830,882. Diagonal scaling and reductions still need to follow the
schedules already specified above.

Reproduce with:

```bash
python3 tools/check_contraction_row_bounds.py --vertices 14 --edges 23 --prime 32783
```

The complete certificate is saved in
`output/accumulator_bounds/valence_refined_V14_E23.json`. No multiplication kernel
or running process is modified by this refinement.
