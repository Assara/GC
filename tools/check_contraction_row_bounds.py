#!/usr/bin/env python3
"""Certify contraction row sums for fixed source V,E using endpoint valences.

Enumerates target degree multisets, not graphs. Non-graphical degree multisets
are retained, so the result is a conservative upper bound.
"""
import argparse
from collections import Counter
from math import comb
import json


def excess_partitions(total, cap, length):
    if total == 0:
        yield []
    elif length:
        for value in range(min(total, cap), 0, -1):
            for tail in excess_partitions(total - value, value, length - 1):
                yield [value] + tail


def row_bound(vertices, edges):
    n = vertices - 1
    excess = 2 * (edges - 1) - 3 * n
    best = {'bound': 0, 'target_valence_counts': {}, 'contributions': [], 'profiles': 0}
    if n < 4 or excess < 0 or 2 * (edges - 1) > n * (n - 1):
        return best
    for partition in excess_partitions(excess, n - 4, n):
        h = Counter([3 + x for x in partition] + [3] * (n - len(partition)))
        total = 0
        contributions = []
        for degree, vertex_count in sorted(h.items()):
            for k in range(2, degree // 2 + 1):
                a, b = k + 1, degree - k + 1
                na, nb = h[a] + 1, h[b] + 1
                if a == b:
                    na = h[a] + 2
                    multiplicity = min(edges, comb(na, 2), a * na // 2)
                    partitions = comb(degree, k) // 2
                else:
                    multiplicity = min(edges, na * nb, a * na, b * nb)
                    partitions = comb(degree, k)
                contribution = vertex_count * partitions * multiplicity
                total += contribution
                contributions.append({
                    'target_valence': degree, 'vertices': vertex_count,
                    'split_endpoint_valences': [a, b],
                    'partitions_per_vertex': partitions,
                    'multiplicity_bound': multiplicity, 'contribution': contribution,
                })
        best['profiles'] += 1
        if total > best['bound']:
            best.update(bound=total, target_valence_counts=dict(sorted(h.items())),
                        contributions=contributions)
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vertices', type=int, default=14)
    parser.add_argument('--edges', type=int, default=23)
    parser.add_argument('--prime', type=int, default=32783)
    args = parser.parse_args()
    v, e, p = args.vertices, args.edges, args.prime
    if v < 4 or p < 2 or not 3*v <= 2*e <= v*(v-1):
        parser.error('require V>=4, p>=2, and 3V<=2E<=V(V-1)')
    down = row_bound(v, e)
    up = row_bound(v + 1, e + 1)
    a = max(down['bound'], e + 1)
    at = e + up['bound']
    print(json.dumps({
        'source_vertices': v, 'source_edges': e, 'prime': p,
        'down_row': down, 'up_row': up,
        'stack_A_row_bound': a, 'stack_AT_row_bound': at,
        'separate_product_accumulator_bound': (p-1)*max(a, at),
        'scope': 'simple graphs, minimum valence three, integer contraction coefficients',
    }, indent=2))


if __name__ == '__main__':
    main()
