#!/usr/bin/env python3
"""Exact integer certificate for GraphHomology/AccumulatorBounds.md.

E is the middle graph's edge count. The contraction window also contains
E-1 and E+1 edge graphs. No floating-point or stored coefficient-width bounds.
"""
import argparse
import json
import math

DEFAULT_PRIME = 32783


def splits_at_vertex(degree):
    return (1 << (degree - 1)) - degree - 1


def split_bound(edges):
    """Convex degree-sequence relaxation for simple, minimum-degree-three graphs."""
    best = 0
    for vertices in range(4, 2 * edges // 3 + 1):
        if 2 * edges > vertices * (vertices - 1):
            continue
        if vertices == 4:  # Only K4 is possible; every degree is three.
            continue
        excess = 2 * edges - 3 * vertices
        full, remainder = divmod(excess, vertices - 4)
        value = full * splits_at_vertex(vertices - 1)
        if remainder:
            value += splits_at_vertex(3 + remainder)
        best = max(best, value)
    return best


def envelope(max_edges):
    values = []
    maximum = 0
    for edges in range(max_edges + 1):
        maximum = max(maximum, split_bound(edges))
        values.append(maximum)
    return values


def bounds(edges, split_envelope, prime):
    q = prime - 1
    r = max(edges + 1, edges * split_envelope[edges - 1])
    t = edges + (edges + 1) * split_envelope[edges]
    return r, t, {
        'each_product': q * max(r, t),
        'each_product_fused_diagonal': q * q * max(r, t),
        'both_products_no_diagonals': q * r * t,
        'both_products_middle_diagonal': q * q * r * t,
        'both_products_both_diagonals': q ** 3 * r * t,
    }


def verify_degree_relaxation():
    # Independent dynamic program over ALL degree lists for small vertex counts.
    # These include non-graphical lists, exactly as the documented relaxation does.
    for vertices in range(4, 11):
        best = {0: 0}
        for _ in range(vertices):
            following = {}
            for degree_sum, score in best.items():
                for degree in range(3, vertices):
                    total = degree_sum + degree
                    following[total] = max(following.get(total, -1),
                                           score + splits_at_vertex(degree))
            best = following
        for total, expected in best.items():
            if vertices == 4:
                actual = 0
            else:
                full, remainder = divmod(total - 3 * vertices, vertices - 4)
                actual = full * splits_at_vertex(vertices - 1)
                if remainder:
                    actual += splits_at_vertex(3 + remainder)
            assert actual == expected, (vertices, total, actual, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edges', type=int, default=23)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--prime', type=int)
    group.add_argument('--prime-bits', type=int, choices=(16, 32),
                       help='certify every prime strictly below 2**BITS')
    args = parser.parse_args()
    if args.prime_bits is not None:
        args.prime = (1 << args.prime_bits) - 1
    elif args.prime is None:
        args.prime = DEFAULT_PRIME
    if args.prime < 2:
        parser.error('--prime must be at least 2 (the arithmetic bound also holds for composite moduli)')
    if not 6 <= args.edges <= 1000:
        parser.error('--edges must be between 6 and 1000')
    verify_degree_relaxation()
    upper = max(256, args.edges)
    split_envelope = envelope(upper)
    r, t, current = bounds(args.edges, split_envelope, args.prime)
    thresholds = {}
    for bits in (32, 64):
        maximum = (1 << (bits - 1)) - 1
        thresholds[bits] = {}
        for scheme in current:
            last = None
            for edges in range(6, upper + 1):
                if bounds(edges, split_envelope, args.prime)[2][scheme] > maximum:
                    break
                last = edges
            else:
                raise RuntimeError('Increase the search range to certify the first failing bound')
            thresholds[bits][scheme] = last
    def integer_root(value, power):
        if power == 1:
            return value
        if power == 2:
            return math.isqrt(value)
        lo, hi = 0, 1 << ((value.bit_length() + power - 1) // power)
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if mid ** power <= value:
                lo = mid
            else:
                hi = mid - 1
        return lo
    forms = {
        'each_product': (max(r, t), 1),
        'each_product_fused_diagonal': (max(r, t), 2),
        'both_products_no_diagonals': (r*t, 1),
        'both_products_middle_diagonal': (r*t, 2),
        'both_products_both_diagonals': (r*t, 3),
    }
    prime_limits = {
        bits: {name: 1 + integer_root(((1 << (bits-1))-1)//factor, power)
               for name, (factor, power) in forms.items()}
        for bits in (32, 64)
    }
    for bits, limits in prime_limits.items():
        maximum = (1 << (bits - 1)) - 1
        for name, prime_limit in limits.items():
            factor, power = forms[name]
            assert (prime_limit - 1) ** power * factor <= maximum
            assert prime_limit ** power * factor > maximum
    print(json.dumps({
        'modulus_bound': args.prime,
        'prime_bits': args.prime_bits,
        'residue_bound': args.prime - 1,
        'middle_edges': args.edges,
        'largest_source_edges': args.edges + 1,
        'split_envelope_lower': split_envelope[args.edges - 1],
        'split_envelope_middle': split_envelope[args.edges],
        'A_absolute_row_sum_bound': r,
        'AT_absolute_row_sum_bound': t,
        'absolute_accumulator_bounds': current,
        'certified_middle_edge_thresholds': thresholds,
        'max_modulus_from_bound_at_middle_edges': prime_limits,
        'signed_limits': {bits: (1 << (bits - 1)) - 1 for bits in (32, 64)},
        'small_degree_relaxation_dp_check': 'passed',
    }, indent=2))


if __name__ == '__main__':
    main()
