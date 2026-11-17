#!/usr/bin/env python3
"""Structural validation of the P7 enumeration and larger-path numerical tools.

Enumeration and Boolean comparisons use exact integers. Gaussian comparisons
are numerical consistency checks and carry no rigorous integration guarantee.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from path_cone import (compact_types, count_types, covariance, enumerate_types,
                       generators, structure, tree_support)
from path_gadgets import cost_table
from path_gradient import cost_gradient, cone_directional_derivatives
from large_path_numerics import (correlation_to_gram_gradient, cost_gradient_qmc,
                                 generator_directions, table_conditional_samples)
from verify_cone import (boundary_labelings, check_finite_covariance, distances,
                         example_graphs, oriented_paths)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def legacy_enumeration(k):
    """Independent original exhaustive rejection algorithm, for small cases."""
    owners, _, _, neighbors = structure(k)
    result = set()
    for increments in itertools.product((-1, 0, 1), repeat=k-1):
        core = [0]
        for delta in increments:
            core.append(core[-1]+delta)
        for offsets in itertools.product((-1, 0, 1), repeat=k+2):
            t = tuple(core)+tuple(core[v]+d for v, d in zip(owners, offsets))
            descending = [min(t[u] for u in neighbors[v]) == t[v]-1 for v in range(k)]
            if all(descending):
                result.add(t)
                continue
            exceptions = [v for v in range(k) if not descending[v]]
            if len(exceptions) != 1:
                continue
            source = exceptions[0]
            low = t[source]
            if any(t[u] != low+1 for u in neighbors[source]):
                continue
            if any(t[v] <= low for v in range(k) if v != source):
                continue
            if any(t[k+j] < low+int(abs(owner-source) <= 1) for j, owner in enumerate(owners)):
                continue
            result.add(t)
    return tuple(sorted(result))


def encode(row):
    value = 0
    for x in row:
        value = (value << 5) | (int(x)+8)
    return value


def verify_p7_rows(types):
    k = 7
    _, edges, _, neighbors = structure(k)
    generic = special = 0
    for start in range(0, len(types), 32768):
        rows = types[start:start+32768].astype(np.int16)
        require(np.all(rows[:, 0] == 0), 'root offsets are not zero')
        for u, v in edges:
            require(np.all(np.abs(rows[:, u]-rows[:, v]) <= 1), 'edge Lipschitz condition failed')
        descending = np.column_stack([np.min(rows[:, neighbors[v]], axis=1) == rows[:, v]-1 for v in range(k)])
        generic += int(np.sum(np.all(descending, axis=1)))
        other = ~np.all(descending, axis=1)
        require(np.all(np.sum(~descending[other], axis=1) == 1), 'special type has multiple exceptions')
        for row, desc in zip(rows[other], descending[other]):
            source = int(np.flatnonzero(~desc)[0])
            low = row[source]
            require(all(row[u] == low+1 for u in neighbors[source]), 'special source neighbor condition failed')
            require(all(row[v] > low for v in range(k) if v != source), 'special source not unique core minimum')
            for j, owner in enumerate(structure(k)[0]):
                require(row[k+j] >= low+int(abs(owner-source) <= 1), 'simple-graph distinctness failed')
            special += 1
    require(generic == count_types(k)['generic'] and special == count_types(k)['special'],
            'enumerated category counts disagree with exact product counts')
    require(not np.any(np.all(types[1:] == types[:-1], axis=1)), 'duplicate types')
    return {'generic_types_checked': generic, 'special_types_checked': special,
            'all_rows_satisfy_original_integer_predicates': True, 'no_duplicate_rows': True}


def verify_p7_graphs(types):
    allowed = {encode(row) for row in types}
    results = []
    for name, adjacency in example_graphs():
        if name not in ('cube', 'Petersen', '7-prism'):
            continue
        dist = distances(adjacency)
        checked = path_count = repeated = chord = 0
        covariance_checks = {}
        for path in oriented_paths(adjacency, 7):
            path_count += 1
            for boundary in boundary_labelings(adjacency, path):
                labels = path+boundary
                rep = len(set(boundary)) < len(boundary)
                ch = any(v in path for v in boundary)
                repeated += int(rep)
                chord += int(ch)
                actual_types = []
                for source in range(len(adjacency)):
                    row = tuple(dist[source][v]-dist[source][path[0]] for v in labels)
                    require(encode(row) in allowed, 'actual P7 source type missing')
                    checked += 1
                    actual_types.append(row)
                category = (rep, ch)
                if category not in covariance_checks:
                    unique = tuple(sorted(set(actual_types)))
                    error, pairs = check_finite_covariance(adjacency, dist, path, boundary,
                                                          unique, generators(unique))
                    covariance_checks[category] = {'max_float_error': error, 'duplicate_edge_pairs': pairs}
        results.append({'graph': name, 'oriented_P7_paths': path_count,
                        'source_patterns_checked': checked,
                        'boundary_labelings_with_repetitions': repeated,
                        'boundary_labelings_with_chords': chord,
                        'independent_covariance_comparisons': list(covariance_checks.values())})
    return results


def numerical_checks():
    support, exact_masses = tree_support(4)
    gen = generators(support)
    masses = np.array([float(x) for x in exact_masses])
    cov = covariance(4, masses, gen)
    ref = cost_gradient(cov, 64)
    normals = ndtri(qmc.Sobol(3, scramble=True, seed=61729).random_base2(17))
    estimate = cost_gradient_qmc(cov, normals)
    gradient_error = float(np.max(np.abs(estimate-ref)))
    require(gradient_error < 2e-4, 'conditional gradient numerical check failed')
    gram = (gen.T*masses) @ gen
    pullback = correlation_to_gram_gradient(4, gram, ref)
    direct = generator_directions(np.asarray(support), gram, pullback)
    reference = cone_directional_derivatives(4, masses, gen, 64)
    adjoint_error = float(np.max(np.abs(direct-reference)))
    require(adjoint_error < 2e-12, 'adjoint normalization calculus disagrees with forward derivative')
    # General table integration must include nonzero tails and mixed signs.
    table = np.full(16, -2, dtype=np.int8)
    values = table_conditional_samples(cov[:4, :4], table, normals)
    constant_error = float(np.max(np.abs(values+2)))
    require(constant_error < 1e-14, 'constant Boolean table integration failed')
    identity = np.eye(4)
    mean = float(np.mean(table_conditional_samples(identity, cost_table(2), normals)))
    require(abs(mean-1.125) < 3e-4, 'independent Boolean expectation check failed')
    return {'P4_conditional_QMC_gradient_max_error_against_deterministic': gradient_error,
            'adjoint_vs_forward_chain_rule_max_error': adjoint_error,
            'constant_table_integration_error': constant_error,
            'independent_P2_table_mean': mean,
            'status': 'numerical diagnostics only, not certified error bounds'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    small = []
    for k in range(2, 6):
        new, old = enumerate_types(k), legacy_enumeration(k)
        require(new == old, 'new enumeration differs from original exhaustive algorithm')
        compact = compact_types(k)
        require(np.array_equal(compact, np.array(old, dtype=np.int8)), 'compact lexicographic order changed')
        small.append({'path_vertices': k, 'types': len(new), 'exact_legacy_equality': True})
    types = compact_types(7)
    result = {'schema_version': 1,
              'status': 'EXACT STRUCTURAL ENUMERATION CHECKS AND NUMERICAL CONSISTENCY CHECKS; NO CONE UPPER BOUND',
              'legacy_comparisons': small, 'P7_counts': count_types(7),
              'P7_compact_bytes': types.nbytes,
              'P7_sorted_int8_sha256': hashlib.sha256(types.tobytes()).hexdigest(),
              'P7_row_validation': verify_p7_rows(types),
              'P7_actual_graph_validation': verify_p7_graphs(types),
              'numerical_validation': numerical_checks(), 'all_checks_passed': True}
    text = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
