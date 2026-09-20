#!/usr/bin/env python3
"""Absolute two-step path bounds for the transposed natural GC Laplacian.

Counts before cancellations; assumes simple graphs, minimum valence three,
no random diagonals inside the composition, and canonical residues modulo p.
"""
import argparse
from functools import lru_cache
import json


@lru_cache(None)
def partitions(total, cap, length):
    if total == 0:
        return ((),)
    if length <= 0 or cap <= 0 or total > cap * length:
        return ()
    return tuple((value,) + tail
                 for value in range(min(total, cap), 0, -1)
                 for tail in partitions(total-value, value, length-1))


def column_bound(vertices, edges):
    excess = 2*edges - 3*vertices
    if vertices < 4 or excess < 0 or 2*edges > vertices*(vertices-1):
        return None
    best = {'bound': -1}
    for profile in partitions(excess, vertices-4, vertices):
        u = [2**t - 1 for t in profile]
        splits = 4*sum(u) - excess
        contract_split = (edges*(splits+3)
                          + 4*sum((t+3)*value for t, value in zip(profile, u))
                          + 4*(sum(u)**2 - sum(value*value for value in u)))
        split_contract = (edges+1)*splits
        bound = contract_split + split_contract
        if bound > best['bound']:
            best = {'bound': bound, 'nontrivalent_valences': [t+3 for t in profile],
                    'trivalent_vertices': vertices-len(profile),
                    'contract_then_split': contract_split,
                    'split_then_contract': split_contract}
    lower_n = vertices-1
    lower_t = 2*(edges-1)-3*lower_n
    lower_splits = max((sum(2**(x+2)-x-4 for x in ps)
                        for ps in partitions(lower_t, lower_n-4, lower_n)), default=0) if lower_t>=0 else 0
    first_stage = max(edges+1, lower_splits)
    return dict(vertices=vertices, edges=edges, first_stage_bound=first_stage, **best)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=32783)
    parser.add_argument('--max-loop', type=int, default=14)
    args = parser.parse_args()
    if args.prime < 2 or not 3 <= args.max_loop <= 20:
        parser.error('require p>=2 and 3<=max-loop<=20')
    q = args.prime-1
    limit = (2**31-1)//q
    loops = []
    for loop in range(3, args.max_loop+1):
        choices = [column_bound(v, loop+v-1) for v in range(4, 2*loop-1)]
        worst = max((x for x in choices if x is not None), key=lambda x: x['bound'])
        worst = dict(loop=loop, **worst,
                     accumulator_bound=q*worst['bound'], certified_int32=worst['bound']<=limit)
        if worst['certified_int32']:
            assert all(x is None or q*x['first_stage_bound'] <= 2**31-1 for x in choices)
        loops.append(worst)
    degree_three = []
    for loop in range(6, args.max_loop+1):
        case = column_bound(loop+4, 2*loop+3)
        if case['bound']<=limit:
            assert q*case['first_stage_bound']<=2**31-1
        degree_three.append(dict(loop=loop, **case,
                                 accumulator_bound=q*case['bound'], certified_int32=case['bound']<=limit))
    assert column_bound(14, 23)['bound'] == 3121
    print(json.dumps({'prime': args.prime, 'max_int32_path_weight': limit,
                      'all_vertex_counts_by_loop': loops,
                      'degree_three_targets': degree_three}, indent=2))


if __name__ == '__main__':
    main()
