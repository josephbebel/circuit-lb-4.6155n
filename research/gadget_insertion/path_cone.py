#!/usr/bin/env python3
"""Finite distance-type relaxation for simultaneous path insertion.

The enumeration is exact integer combinatorics.  Covariance calculations use
ordinary floating point and are exploration tools, NOT interval certificates.
The cone contains every actual path neighborhood (up to a vanishing truncation
shell); it also contains patterns not necessarily realized by any cubic graph.
See PROOF.md for the remaining universal certification obligations.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json

import numpy as np


def boundary_owners(k: int) -> tuple[int, ...]:
    if not isinstance(k, int) or k < 2:
        raise ValueError("a path must have at least two vertices")
    return (0, 0, *range(1, k - 1), k - 1, k - 1)


def structure(k: int):
    owners = boundary_owners(k)
    edges = [(v, v + 1) for v in range(k - 1)]
    boundary = [(owner, k + j) for j, owner in enumerate(owners)]
    edges += boundary
    neighbors = [[] for _ in range(k)]
    for u, v in edges:
        if u < k:
            neighbors[u].append(v)
        if v < k:
            neighbors[v].append(u)
    assert all(len(row) == 3 for row in neighbors)
    return owners, edges, boundary, neighbors


def enumerate_types(k: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate a complete, possibly overinclusive set of relative distances.

    Core vertices are 0,...,k-1; leaf k+j is attached to owners[j].
    Labels denote occurrences: external neighbors can coincide or be other core
    vertices.  Ignoring such identifications enlarges the cone safely.

    A generic source lies outside the core.  A special source is one core vertex;
    simple-graph constraints exclude a boundary edge duplicating a path edge.
    """
    if k > 6:
        raise ValueError("explicit enumeration is limited to k <= 6")
    owners, _, _, neighbors = structure(k)
    result = set()
    leaf_choices = tuple(itertools.product((-1, 0, 1), repeat=k + 2))
    for increments in itertools.product((-1, 0, 1), repeat=k - 1):
        core = [0]
        for delta in increments:
            core.append(core[-1] + delta)
        for leaf_delta in leaf_choices:
            t = tuple(core) + tuple(core[v] + d for v, d in zip(owners, leaf_delta))
            descends = [min(t[u] for u in neighbors[v]) == t[v] - 1
                        for v in range(k)]
            if all(descends):
                result.add(t)
                continue
            # If z is a core vertex, that vertex is the only one whose shortest
            # path to z does not start with a step down in distance.
            exceptions = [v for v in range(k) if not descends[v]]
            if len(exceptions) != 1:
                continue
            source = exceptions[0]
            minimum = t[source]
            if any(t[u] != minimum + 1 for u in neighbors[source]):
                continue
            if any(t[v] <= minimum for v in range(k) if v != source):
                continue
            good = True
            for j, owner in enumerate(owners):
                # A boundary neighbor cannot be its owner or an adjacent core
                # vertex (that would be a repeated required path edge).
                distinct = abs(owner - source) <= 1
                if t[k + j] < minimum + int(distinct):
                    good = False
                    break
            if good:
                result.add(t)
    return tuple(sorted(result))


def extended_tree_distances(k: int) -> np.ndarray:
    owners = boundary_owners(k)
    n = 2 * k + 2
    distances = np.zeros((n, n), dtype=np.int64)
    for x in range(n):
        owner_x, extra_x = (x, 0) if x < k else (owners[x - k], 1)
        for y in range(x):
            owner_y, extra_y = (y, 0) if y < k else (owners[y - k], 1)
            d = abs(owner_x - owner_y) + extra_x + extra_y
            distances[x, y] = distances[y, x] = d
    return distances


def tree_distribution(k: int, types=None):
    """Exact rational masses of the k+2 outward branches of the infinite tree."""
    types = enumerate_types(k) if types is None else types
    lookup = {t: i for i, t in enumerate(types)}
    distances = extended_tree_distances(k)
    result = [Fraction(0) for _ in types]
    for j, owner in enumerate(boundary_owners(k)):
        leaf = k + j
        t = tuple(int(x) for x in distances[:, leaf] - distances[0, leaf])
        result[lookup[t]] += Fraction(1, 3 * 2**owner)
    assert sum(result) == 1
    return result


def generators(types) -> np.ndarray:
    return np.exp2(-np.asarray(types, dtype=float) / 2)


def correlation(gram: np.ndarray) -> np.ndarray:
    diagonal = np.diag(gram)
    if np.any(diagonal <= 0):
        raise ValueError("positive diagonal required")
    roots = np.sqrt(diagonal)
    return gram / np.outer(roots, roots)


def covariance(k: int, masses, gen: np.ndarray) -> np.ndarray:
    """Ordinary floating approximation of corr(B corr(G(m)) B^T)."""
    masses = np.asarray(masses, dtype=float)
    if masses.shape != (len(gen),) or not np.all(np.isfinite(masses)):
        raise ValueError("invalid mass vector")
    if np.any(masses < 0) or masses.sum() <= 0:
        raise ValueError("nonnegative, nonzero mass vector required")
    masses = masses / masses.sum()
    gram = (gen.T * masses) @ gen
    vertex_correlation = correlation(gram)
    owners = boundary_owners(k)
    incidence = np.zeros((k + 2, 2 * k + 2))
    for j, owner in enumerate(owners):
        incidence[j, owner] = incidence[j, k + j] = 1
    return correlation(incidence @ vertex_correlation @ incidence.T)


def check(k: int):
    types = enumerate_types(k)
    gen = generators(types)
    masses = tree_distribution(k, types)
    cov = covariance(k, masses, gen)
    owners = boundary_owners(k)
    r = lambda d: (1 + d / 3) * 2**(-d / 2)
    expected = np.eye(k + 2)
    for i, u in enumerate(owners):
        for j, v in enumerate(owners[:i]):
            d = abs(u - v)
            expected[i, j] = expected[j, i] = (r(d) + 2*r(d+1) + r(d+2)) / (2+2*r(1))
    error = float(np.max(np.abs(cov - expected)))
    assert error < 5e-14, error
    if k == 2:
        assert len(types) == 123
    return {
        "path_vertices": k,
        "boundary_occurrences": k + 2,
        "distance_types": len(types),
        "generic_types": sum(all(min(t[u] for u in structure(k)[3][v]) == t[v]-1
                                 for v in range(k)) for t in types),
        "tree_support_size": sum(bool(x) for x in masses),
        "tree_mass_sum": str(sum(masses)),
        "tree_covariance_max_float_error": error,
        "status": "integer enumeration and numerical consistency check; no cone maximum certified",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[2, 3, 4])
    args = parser.parse_args()
    for size in args.sizes:
        print(json.dumps(check(size), sort_keys=True), flush=True)
