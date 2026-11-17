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
import math

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


def _boundary_option_blocks(k: int):
    """Yield disjoint (core, boundary-option blocks, category) products.

    A core vertex with an existing downward path neighbor needs no further
    descent constraint on its boundary leaves. Otherwise at least one such
    leaf must be one level lower. For a core source, that source is the unique
    minimum among core levels, its neighbors are one level higher, and the
    simple-graph distinctness restrictions are imposed directly on leaf choices.
    This factors the old exhaustive rejection loop without weakening any test.
    """
    boundary_owners(k)
    if k > 7:
        raise ValueError("explicit enumeration is limited to k <= 7")
    for increments in itertools.product((-1, 0, 1), repeat=k - 1):
        core = [0]
        for delta in increments:
            core.append(core[-1] + delta)
        path_neighbors = [(core[v-1:v] if v else []) + core[v+1:v+2]
                          for v in range(k)]

        def options(source):
            blocks = []
            for v in range(k):
                count = 2 if v in (0, k-1) else 1
                if v == source:
                    choices = ((core[v] + 1,) * count,)
                else:
                    choices = tuple(itertools.product(range(core[v]-1, core[v]+2), repeat=count))
                    if min(path_neighbors[v]) != core[v]-1:
                        choices = tuple(row for row in choices if min(row) == core[v]-1)
                    if source is not None:
                        floor = core[source] + int(abs(v-source) <= 1)
                        choices = tuple(row for row in choices if min(row) >= floor)
                if not choices:
                    return None
                blocks.append(choices)
            return tuple(blocks)

        generic = options(None)
        if generic is not None:
            yield tuple(core), generic, "generic"
        minimum = min(core)
        if core.count(minimum) == 1:
            source = core.index(minimum)
            if all(value == minimum+1 for value in path_neighbors[source]):
                special = options(source)
                if special is not None:
                    yield tuple(core), special, "special"


def iter_types(k: int):
    """Stream each type exactly once; order is deterministic but not sorted."""
    for core, blocks, _ in _boundary_option_blocks(k):
        for leaves in itertools.product(*blocks):
            yield core + tuple(itertools.chain.from_iterable(leaves))


def count_types(k: int) -> dict[str, int]:
    """Count the disjoint Cartesian products using exact integer arithmetic."""
    result = {"generic": 0, "special": 0}
    for _, blocks, category in _boundary_option_blocks(k):
        result[category] += math.prod(len(block) for block in blocks)
    result["total"] = result["generic"] + result["special"]
    return result


def compact_types(k: int, sort=True) -> np.ndarray:
    """Materialize types as int8 rows (about 16.6 MB for all P7 types).

    Sorted output has the same lexicographic indices as enumerate_types.
    Set sort=False when only streaming indices or actual type vectors matter.
    """
    expected = count_types(k)["total"]
    result = np.empty((expected, 2*k+2), dtype=np.int8)
    count = 0
    for row in iter_types(k):
        if count >= expected:
            raise ArithmeticError("type iterator exceeded its exact product count")
        result[count] = row
        count += 1
    if count != expected:
        raise ArithmeticError("type iterator did not reach its exact product count")
    if sort:
        order = np.lexsort(tuple(result[:, col] for col in range(2*k+1, -1, -1)))
        result = result[order]
    return result


def enumerate_types(k: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate a complete, possibly overinclusive set of relative distances.

    Core vertices are 0,...,k-1; leaf k+j is attached to owners[j].
    Labels denote occurrences: external neighbors can coincide or be other core
    vertices.  Ignoring such identifications enlarges the cone safely.

    A generic source lies outside the core.  A special source is one core vertex;
    simple-graph constraints exclude a boundary edge duplicating a path edge.

    This compatibility API retains the original sorted tuple representation.
    For large cones prefer compact_types or iter_types to avoid tuple overhead.
    """
    return tuple(sorted(iter_types(k)))


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
    lookup = {tuple(t): i for i, t in enumerate(types)}
    result = [Fraction(0) for _ in types]
    support, weights = tree_support(k)
    for t, weight in zip(support, weights):
        result[lookup[t]] += weight
    assert sum(result) == 1
    return result


def tree_support(k: int):
    """Return the k+2 tree types and exact masses, without enumerating the cone.

    Rows are ordered by their boundary branch, not lexicographically.
    """
    distances = extended_tree_distances(k)
    support, weights = [], []
    for j, owner in enumerate(boundary_owners(k)):
        leaf = k+j
        support.append(tuple(int(x) for x in distances[:, leaf]-distances[0, leaf]))
        weights.append(Fraction(1, 3*2**owner))
    if sum(weights) != 1:
        raise ArithmeticError("tree masses do not sum to one")
    return tuple(support), tuple(weights)


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
    types = compact_types(k)
    support, masses = tree_support(k)
    gen = generators(support)
    for row in support:
        if not np.any(np.all(types == np.asarray(row), axis=1)):
            raise ArithmeticError("a tree support type is absent")
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
        "generic_types": count_types(k)["generic"],
        "special_types": count_types(k)["special"],
        "compact_type_storage_bytes": types.nbytes,
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
