#!/usr/bin/env python3
"""Exhaustive finite checks of the Boolean path rules, not a Gaussian proof."""
import argparse
import json

import numpy as np

from path_gadgets import (ConditionalGaussianCost, Qsqrt2, common_tree_covariance,
                          cost_table, exact_tree_covariance, least_minimizer,
                          prefix_indices, tree_covariance)


def objective(boundary, internal):
    """Direct definition, independent of the optimizer's dynamic program."""
    n = len(internal) + 1
    count = int(len({int(boundary[0]), int(boundary[1]), int(internal[0])}) > 1)
    for v in range(1, n - 1):
        count += int(len({int(internal[v - 1]), int(internal[v]), int(boundary[v + 1])}) > 1)
    return count + int(len({int(internal[-1]), int(boundary[-2]), int(boundary[-1])}) > 1)


def verify(path_size):
    b = path_size + 2
    masks = np.arange(1 << b, dtype=np.uint64)
    bits = ((masks[:, None] >> np.arange(b, dtype=np.uint64)) & 1).astype(bool)
    internal_masks = np.arange(1 << (path_size - 1), dtype=np.uint64)
    internal = ((internal_masks[:, None] >> np.arange(path_size - 1, dtype=np.uint64)) & 1).astype(bool)
    table = cost_table(path_size)
    least = least_minimizer(bits)
    for mask, boundary in enumerate(bits):
        direct = np.array([objective(boundary, x) for x in internal])
        best = direct.min()
        assert table[mask] == best, ("DP mismatch", path_size, mask)
        minimizers = internal[direct == best]
        intersection = np.logical_and.reduce(minimizers)
        assert objective(boundary, intersection) == best, ("lattice least not optimal", path_size, mask)
        assert np.array_equal(intersection, least[mask]), ("least-minimizer mismatch", path_size, mask)
    assert np.array_equal(table, table[::-1]), ("complement symmetry", path_size)
    submodular_checks = monotone_checks = 0
    for i in range(b):
        base = masks[(masks & (1 << i)) == 0].astype(np.int64)
        assert np.all(least[base] <= least[base | (1 << i)]), ("monotonicity", path_size, i)
        monotone_checks += len(base)
        for j in range(i):
            base_ij = base[(base & (1 << j)) == 0]
            mixed = table[base_ij | (1 << i) | (1 << j)] - table[base_ij | (1 << i)] - table[base_ij | (1 << j)] + table[base_ij]
            assert np.all(mixed <= 0), ("submodularity", path_size, i, j)
            submodular_checks += len(base_ij)
    return {"path_vertices": path_size, "boundary_patterns": 1 << b,
            "internal_assignments_per_pattern": 1 << (path_size - 1),
            "monotonicity_checks": monotone_checks, "submodularity_checks": submodular_checks,
            "boolean_checks_passed": True}


def check_covariances(max_size=12):
    from fractions import Fraction
    c2 = exact_tree_covariance(2)
    assert c2[0][1] == Qsqrt2(Fraction(1, 4), Fraction(1, 2))
    assert c2[0][2] == Qsqrt2(Fraction(1, 2), Fraction(1, 4))
    common = common_tree_covariance(max_size)
    for n in range(2, max_size + 1):
        exact = exact_tree_covariance(n)
        assert all(exact[i][i] == Qsqrt2(Fraction(1)) for i in range(n + 2))
        assert all(exact[i][j] == exact[j][i] for i in range(n + 2) for j in range(n + 2))
        idx = prefix_indices(n)
        assert np.allclose(common[np.ix_(idx, idx)], tree_covariance(n), atol=2e-14, rtol=0)
    one = ConditionalGaussianCost(np.ones((4, 4)))
    # The floating eigensolver can introduce residual factors around 1e-8.
    assert abs(one.sample_values(np.zeros((1, 3)))[0]) < 1e-14
    return {"covariance_checks_passed": True, "max_path_vertices": max_size,
            "note": "Algebraic entries checked exactly; matrix factorization remains floating point."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-size', type=int, default=7)
    args = parser.parse_args()
    if not 2 <= args.max_size <= 9:
        parser.error('--max-size must be in [2, 9]')
    for n in range(2, args.max_size + 1):
        print(json.dumps(verify(n)), flush=True)
    print(json.dumps(check_covariances()), flush=True)


if __name__ == '__main__':
    main()
