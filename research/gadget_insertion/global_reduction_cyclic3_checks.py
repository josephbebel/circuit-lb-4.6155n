#!/usr/bin/env python3
"""Exact finite checks for GLOBAL_REDUCTION_CYCLIC3.md."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from global_reduction_checks import Graph, generalized_petersen, random_cubic, require


def cuts3(g):
    """Enumerate minimal 3-edge cuts of a 3-edge-connected graph."""
    g.check_class()
    require(not g.two_cuts(), 'input is not 3-edge-connected')
    rows = g.adjacency()
    all_vertices = frozenset(range(g.n))
    found = {}
    for cut in itertools.combinations(range(len(g.edges)), 3):
        omitted = set(cut)
        seen = {0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v, edge in rows[u]:
                if edge not in omitted and v not in seen:
                    seen.add(v)
                    stack.append(v)
        if len(seen) == g.n:
            continue
        a = frozenset(seen)
        b = all_vertices - a
        require(all((g.edges[e][0] in a) != (g.edges[e][1] in a) for e in cut),
                'enumerated cut is not minimal')
        a = min((a, b), key=lambda x: (len(x), tuple(sorted(x))))
        found[a] = cut
    return found


def complete(g, shore):
    vertices = sorted(shore)
    local = {v: i for i, v in enumerate(vertices)}
    edges = []
    for u, v in g.edges:
        if u in shore and v in shore:
            edges.append((local[u], local[v]))
        elif u in shore:
            edges.append((local[u], len(vertices)))
        elif v in shore:
            edges.append((local[v], len(vertices)))
    h = Graph(len(vertices) + 1, tuple(edges))
    h.check_class()
    require(not h.two_cuts(), 'completion has a 2-edge cut')
    return h, vertices


def sum3(g, h, u=0, v=0):
    gu = [x for x, _ in g.adjacency()[u]]
    hv = [x for x, _ in h.adjacency()[v]]
    labels = [(0, x) for x in range(g.n) if x != u]
    labels += [(1, x) for x in range(h.n) if x != v]
    local = {x: i for i, x in enumerate(labels)}
    edges = [(local[0, x], local[0, y]) for x, y in g.edges if u not in (x, y)]
    edges += [(local[1, x], local[1, y]) for x, y in h.edges if v not in (x, y)]
    edges += [(local[0, x], local[1, y]) for x, y in zip(gu, hv)]
    answer = Graph(len(labels), tuple(edges))
    answer.check_class()
    require(not answer.two_cuts(), '3-sum lost 3-edge-connectivity')
    return answer


def replace_all_vertices(g, h, remove=0):
    ports = [v for v, _ in h.adjacency()[remove]]
    require(len(set(ports)) == 3, 'replacement needs distinct ports')
    fragment = [v for v in range(h.n) if v != remove]
    local = {(a, v): a * len(fragment) + i for a in range(g.n)
             for i, v in enumerate(fragment)}
    edges = [(local[a, u], local[a, v]) for a in range(g.n)
             for u, v in h.edges if remove not in (u, v)]
    used = [0] * g.n
    for a, b in g.edges:
        edges.append((local[a, ports[used[a]]], local[b, ports[used[b]]]))
        used[a] += 1
        used[b] += 1
    answer = Graph(g.n * len(fragment), tuple(edges))
    answer.check_class()
    require(not answer.two_cuts(), 'replacement lost 3-edge-connectivity')
    return answer


def small_core(g, cut_list, threshold):
    n = g.n
    require(n > 3 * threshold, 'outside size regime')
    require(all(min(len(a), n - len(a)) < threshold for a in cut_list),
            'balanced cut present')
    vertices = frozenset(range(n))
    small = {s for a in cut_list for s in (a, vertices - a) if len(s) < threshold}
    parts = sorted((a for a in small if not any(a < b for b in small)),
                   key=lambda a: tuple(sorted(a)))
    require(set().union(*parts) == set(vertices), 'parts do not cover')
    require(sum(map(len, parts)) == n, 'parts overlap')
    owner = {v: i for i, a in enumerate(parts) for v in a}
    edges = [(owner[u], owner[v]) for u, v in g.edges if owner[u] != owner[v]]
    core = Graph(len(parts), tuple(edges))
    core_cuts = cuts3(core)
    require(all(min(len(a), core.n - len(a)) == 1 for a in core_cuts),
            'core has a nontrivial 3-edge cut')
    layouts = []
    for order in (list(range(core.n)), list(reversed(range(core.n))), core.bfs_order()):
        expanded = [v for i in order for v in sorted(parts[i])]
        w0, w1 = core.width(order), g.width(expanded)
        require(2 * (w1 - w0) <= 3 * threshold, 'expansion width bound failed')
        layouts.append({'core_width': w0, 'expanded_width': w1})
    return {'T': threshold, 'core_order': core.n,
            'largest_part': max(map(len, parts)), 'layouts': layouts}


def run():
    k4 = Graph(4, tuple(itertools.combinations(range(4), 2)))
    cube = generalized_petersen(4, 1)
    petersen = generalized_petersen(5, 2)
    families = [('K4', k4), ('cube', cube), ('Petersen', petersen)]
    families += [(f'prism_{n}', generalized_petersen(n, 1)) for n in range(3, 9)]
    for n in range(8, 25, 2):
        g = random_cubic(n, 431 + n)
        if not g.two_cuts():
            families.append((f'random_{n}', g))
    families += [('sum_K4_cube', sum3(k4, cube)),
                 ('sum_Petersen_cube', sum3(petersen, cube)),
                 ('sum_three_cubes', sum3(sum3(cube, cube), cube)),
                 ('cube_triangle_replacements', replace_all_vertices(cube, k4)),
                 ('cube_cube_replacements', replace_all_vertices(cube, cube))]
    results = []
    completion_count = 0
    core_count = 0
    for name, g in families:
        cuts = cuts3(g)
        all_vertices = frozenset(range(g.n))
        for a in cuts:
            b = all_vertices - a
            ha, va = complete(g, a)
            hb, vb = complete(g, b)
            oa, ob = ha.bfs_order(), hb.bfs_order()
            expanded = [va[x] for x in oa if x < len(va)]
            expanded += [vb[x] for x in ob if x < len(vb)]
            require(g.width(expanded) <= max(ha.width(oa), hb.width(ob)) + 3,
                    'completed-shore concatenation failed')
            completion_count += 2
        cores = []
        for threshold in range(2, (g.n - 1) // 3 + 1):
            if all(min(len(a), g.n - len(a)) < threshold for a in cuts):
                cores.append(small_core(g, cuts, threshold))
        core_count += len(cores)
        results.append({'name': name, 'order': g.n, 'three_edge_cuts': len(cuts),
                        'nontrivial_three_edge_cuts': sum(len(a) > 1 for a in cuts),
                        'small_shore_checks': cores})
    return {'status': 'passed', 'graph_cases': len(results),
            'completed_shores': completion_count, 'small_shore_core_cases': core_count,
            'scope': 'Exact finite construction checks; universal proof is in the note.',
            'results': results}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'results'}, indent=2))
