#!/usr/bin/env python3
"""Independent structural checks for the exploratory path distance cones.

Exact checks: comparison with upstream's 123 types; every source vertex for
every oriented simple P2/P3/P4 in several cubic graphs; tree masses and vertex
Gram matrices in Q(sqrt(2)).  Floating covariance comparisons are consistency
checks only.  None of these checks certifies a maximum over the entire cone.
"""

from __future__ import annotations

import argparse
import ast
from collections import deque
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

import path_cone as target


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def upstream_type_set():
    """Execute only the original pure-integer enumerator, avoiding SciPy imports."""
    root = Path(__file__).resolve().parents[2]
    path = root / "mtree" / "enclosure" / "gc4_lib.py"
    data = path.read_bytes()
    digest = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    require(digest == "aef45c027866a30eb5b4312dc738deb939c44054",
            "upstream enumerator differs from the reviewed source version")
    module = ast.parse(data.decode("utf-8"))
    nodes = [node for node in module.body
             if isinstance(node, ast.FunctionDef) and node.name == "enum_types"]
    require(len(nodes) == 1, "upstream enum_types function not unique")
    namespace = {"itertools": itertools,
                 "EDG": [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]}
    code = compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec")
    exec(code, namespace)
    result = namespace["enum_types"]()
    require(len(result) == len({row for _, row in result}) == 123,
            "upstream did not produce exactly 123 distinct types")
    return {row for _, row in result}


def graph(n, edges):
    adjacency = [set() for _ in range(n)]
    for u, v in edges:
        require(0 <= u < n and 0 <= v < n and u != v, "invalid graph edge")
        adjacency[u].add(v)
        adjacency[v].add(u)
    require(all(len(row) == 3 for row in adjacency), "graph is not cubic")
    return tuple(tuple(sorted(row)) for row in adjacency)


def prism(n):
    require(n >= 3, "prism needs at least three vertices on each cycle")
    edges = [(i, (i+1) % n) for i in range(n)]
    edges += [(n+i, n+(i+1) % n) for i in range(n)]
    edges += [(i, n+i) for i in range(n)]
    return graph(2*n, edges)


def example_graphs():
    yield "K4", graph(4, itertools.combinations(range(4), 2))
    yield "K3,3", graph(6, itertools.product(range(3), range(3, 6)))
    yield "cube", graph(8, ((v, v ^ (1 << bit)) for v in range(8) for bit in range(3)))
    petersen = [(i, (i+1) % 5) for i in range(5)]
    petersen += [(5+i, 5+(i+2) % 5) for i in range(5)]
    petersen += [(i, i+5) for i in range(5)]
    yield "Petersen", graph(10, petersen)
    for n in (3, 5, 6, 7):
        yield f"{n}-prism", prism(n)


def distances(adjacency):
    result = []
    for source in range(len(adjacency)):
        row = [-1] * len(adjacency)
        row[source] = 0
        todo = deque([source])
        while todo:
            v = todo.popleft()
            for w in adjacency[v]:
                if row[w] < 0:
                    row[w] = row[v]+1
                    todo.append(w)
        require(min(row) >= 0, "graph is disconnected")
        result.append(row)
    return result


def oriented_paths(adjacency, k):
    def extend(path):
        if len(path) == k:
            yield path
            return
        for v in adjacency[path[-1]]:
            if v not in path:
                yield from extend(path+(v,))
    for v in range(len(adjacency)):
        yield from extend((v,))


def boundary_labelings(adjacency, path):
    """Remove only chosen path edges: retain chords and repeated neighbors."""
    choices = []
    for i, v in enumerate(path):
        removed = set()
        if i:
            removed.add(path[i-1])
        if i+1 < len(path):
            removed.add(path[i+1])
        boundary = tuple(w for w in adjacency[v] if w not in removed)
        require(len(boundary) == (2 if i in (0, len(path)-1) else 1),
                "incorrect number of boundary occurrences")
        choices.append(tuple(itertools.permutations(boundary)))
    for blocks in itertools.product(*choices):
        yield tuple(v for block in blocks for v in block)


def check_finite_covariance(adjacency, dist, path, boundary, types, gen):
    """Compare an actual graph Gaussian field to its distance-cone mixture."""
    n, k = len(adjacency), len(path)
    labels = path+boundary
    lookup = {t: j for j, t in enumerate(types)}
    masses = np.zeros(len(types))
    for source in range(n):
        t = tuple(dist[source][v]-dist[source][path[0]] for v in labels)
        masses[lookup[t]] += 2.0**(-dist[source][path[0]])
    cone = target.covariance(k, masses, gen)

    # Independent construction on actual graph vertices, with repeated labels
    # represented by the same row of the underlying field.
    field = np.exp2(-np.asarray(dist, dtype=float)/2)
    field /= np.linalg.norm(field, axis=1)[:, None]
    owners = (0, 0, *range(1, k-1), k-1, k-1)
    keys = np.asarray([field[path[owner]]+field[v]
                       for owner, v in zip(owners, boundary)])
    keys /= np.linalg.norm(keys, axis=1)[:, None]
    actual = keys @ keys.T
    error = float(np.max(np.abs(actual-cone)))
    require(error < 2e-13, f"finite graph covariance mismatch: {error}")

    same_edge_pairs = 0
    physical_edges = [frozenset((path[owner], v)) for owner, v in zip(owners, boundary)]
    for i, edge in enumerate(physical_edges):
        for j in range(i):
            if edge == physical_edges[j]:
                same_edge_pairs += 1
                require(abs(cone[i, j]-1) < 2e-13,
                        "two occurrences of a chord do not have identical keys")
    return error, same_edge_pairs


def check_actual_graph(name, adjacency, k, types):
    allowed = set(types)
    gen = target.generators(types)
    dist = distances(adjacency)
    paths = labelings = sources = special = repeated = chords = 0
    sampled = set()
    covariance_error = 0.0
    same_edge_pairs = 0
    for path in oriented_paths(adjacency, k):
        paths += 1
        for boundary in boundary_labelings(adjacency, path):
            labelings += 1
            labels = path+boundary
            has_repetition = len(set(boundary)) < len(boundary)
            has_chord = any(v in path for v in boundary)
            repeated += int(has_repetition)
            chords += int(has_chord)
            for source in range(len(adjacency)):
                t = tuple(dist[source][v]-dist[source][path[0]] for v in labels)
                require(t in allowed,
                        f"missing type: graph={name}, path={path}, "
                        f"boundary={boundary}, source={source}, type={t}")
                sources += 1
                special += int(source in path)
            category = (has_repetition, has_chord)
            if category not in sampled:
                error, pairs = check_finite_covariance(
                    adjacency, dist, path, boundary, types, gen)
                covariance_error = max(covariance_error, error)
                same_edge_pairs += pairs
                sampled.add(category)
    return {
        "graph": name, "vertices": len(adjacency), "path_vertices": k,
        "oriented_paths": paths, "boundary_labelings": labelings,
        "source_distance_types_checked": sources,
        "core_sources_checked": special,
        "labelings_with_repeated_boundary_vertices": repeated,
        "labelings_with_chord_occurrences": chords,
        "sampled_covariance_patterns": len(sampled),
        "sampled_duplicate_edge_pairs": same_edge_pairs,
        "covariance_max_float_error": covariance_error,
    }


# The pair (a,b) means a+b*sqrt(2), with exact rational a,b.
def add(x, y):
    return (x[0]+y[0], x[1]+y[1])


def scale(x, q):
    return (x[0]*q, x[1]*q)


def half_power(exponent):
    if exponent % 2 == 0:
        return (Fraction(2)**(exponent//2), Fraction(0))
    return (Fraction(0), Fraction(2)**((exponent-1)//2))


def multiply(x, y):
    return (x[0]*y[0]+2*x[1]*y[1], x[0]*y[1]+x[1]*y[0])


def divide(x, y):
    norm = y[0]*y[0]-2*y[1]*y[1]
    require(norm != 0, "division by zero in Q(sqrt(2))")
    return scale(multiply(x, (y[0], -y[1])), 1/norm)


def exact_float(x):
    return float(x[0])+float(x[1])*np.sqrt(2.0)


def check_tree(k, types):
    owners = (0, 0, *range(1, k-1), k-1, k-1)
    n = 2*k+2
    # Build and traverse the finite extended tree independently of the target's
    # explicit distance formula. Unlike the finite graph tests, boundary leaves
    # here have degree one, so use a direct adjacency construction.
    adjacency = [set() for _ in range(n)]
    edges = [(i, i+1) for i in range(k-1)]+[(owner, k+j) for j, owner in enumerate(owners)]
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    dist = distances(adjacency)
    require(np.array_equal(np.asarray(dist), target.extended_tree_distances(k)),
            "extended-tree distances disagree with BFS")
    masses = target.tree_distribution(k, types)
    require(all(isinstance(w, Fraction) and w >= 0 for w in masses),
            "tree masses are not nonnegative rationals")
    require(sum(masses) == Fraction(1), "tree masses do not sum to exactly one")
    support = [(t, w) for t, w in zip(types, masses) if w]
    require(len(support) == k+2, "unexpected tree support")
    zero = (Fraction(0), Fraction(0))
    gram = [[zero for _ in range(n)] for _ in range(n)]
    for x in range(n):
        for y in range(n):
            value = zero
            for t, weight in support:
                value = add(value, scale(half_power(-t[x]-t[y]), weight))
            expected = scale(half_power(-dist[x][y]), 1+Fraction(dist[x][y], 3))
            require(value == expected,
                    f"exact tree Gram mismatch at k={k}, entry={(x,y)}")
            gram[x][y] = value
    norm = add((Fraction(2), Fraction(0)), scale(gram[0][1], 2))
    exact_keys = np.zeros((k+2, k+2))
    for i, u in enumerate(owners):
        for j, v in enumerate(owners):
            total = zero
            for x in (u, k+i):
                for y in (v, k+j):
                    total = add(total, gram[x][y])
            exact_keys[i, j] = exact_float(divide(total, norm))
    cov = target.covariance(k, masses, target.generators(types))
    error = float(np.max(np.abs(cov-exact_keys)))
    require(error < 2e-13, "tree boundary covariance float mismatch")
    return {"path_vertices": k, "types": len(types),
            "exact_mass_sum": str(sum(masses)), "support_size": len(support),
            "exact_Q_sqrt2_Gram_entries_checked": n*n,
            "boundary_covariance_max_float_error": error}


def verify():
    sizes = (2, 3, 4)
    types = {k: target.enumerate_types(k) for k in sizes}
    require(set(types[2]) == upstream_type_set(),
            "P2 cone does not exactly equal the upstream 123-type cone")
    tree = [check_tree(k, types[k]) for k in sizes]
    examples = [check_actual_graph(name, adjacency, k, types[k])
                for name, adjacency in example_graphs() for k in sizes]
    require(sum(x["labelings_with_repeated_boundary_vertices"] for x in examples) > 0,
            "test set never exercised repeated boundary vertices")
    require(sum(x["sampled_duplicate_edge_pairs"] for x in examples) > 0,
            "test set never exercised the two occurrences of a chord")
    return {
        "checks_passed": True,
        "upstream_123_type_set_exact_match": True,
        "source_distance_types_checked": sum(x["source_distance_types_checked"] for x in examples),
        "tree_checks": tree,
        "finite_graph_checks": examples,
        "scope": "Exact structural tests and rational tree identities; numerical covariance consistency only; no universal cone maximum certified.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(verify(), indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result, encoding="utf-8")
    print(result, end="")


if __name__ == "__main__":
    main()
