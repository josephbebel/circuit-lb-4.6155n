#!/usr/bin/env python3
"""Independent checks of the exploratory Plackett derivative implementation.

Exact Walsh-table inversion is a finite integer check. Everything involving
Gaussian integrals, gradients, finite differences, or stationarity here is a
numerical diagnostic, not an interval enclosure or universal graph bound.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np

from path_cone import covariance, enumerate_types, generators, tree_distribution
from path_gadgets import cost_table, tree_covariance
from path_gradient import (cone_directional_derivatives, cost_gradient,
                           expected_cost_quadrature, fourth_sign_moment,
                           mixed_difference_fourier)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def check_walsh():
    identities = 0
    for k in (2, 3, 4):
        table = cost_table(k)
        for i, j in itertools.combinations(range(k + 2), 2):
            rest, nums, denominator = mixed_difference_fourier(k, i, j)
            for assignment in range(1 << k):
                base = sum(1 << coordinate for pos, coordinate in enumerate(rest)
                           if (assignment >> pos) & 1)
                direct = int(table[base | (1 << i) | (1 << j)]) - int(table[base | (1 << i)]) - int(table[base | (1 << j)]) + int(table[base])
                reconstructed = sum(n * (-1) ** ((mask & assignment).bit_count())
                                    for mask, n in enumerate(nums))
                require(reconstructed == denominator * direct, "Walsh reconstruction failed")
                require(direct <= 0, "a mixed difference was positive")
                identities += 1
    return {"exact_walsh_inversion_identities": identities,
            "all_mixed_differences_nonpositive": True}


def load_upstream():
    root = Path(__file__).resolve().parents[2]
    upstream = root / 'mtree' / 'enclosure'
    sys.path.insert(0, str(upstream))
    import bk7_core
    import gc4_lib
    import phifun
    return bk7_core, gc4_lib, phifun


def check_p2_and_fourth_moments(order=64):
    bk7, gc4, phifun = load_upstream()
    gram_factor = np.array([[3., 1., 0., 1.], [1., 2., -1., 0.],
                            [0., 1., 3., 1.], [1., 0., 1., 2.]])
    gram = gram_factor @ gram_factor.T + np.eye(4)
    scales = np.sqrt(np.diag(gram))
    generic = gram / np.outer(scales, scales)
    tree = tree_covariance(2)
    points = {'tree': tree, 'generic_pd': generic}
    rows = []
    for name, cov in points.items():
        entries = gc4.cvec(cov)
        computed = expected_cost_quadrature(cov, order)
        reference = gc4.phiF(entries)
        require(abs(computed - reference) < 2e-10, 'P2 upstream cost mismatch')
        grad = cost_gradient(cov, order)
        ours = np.array([grad[i, j] for i, j in gc4.PAIRS])
        theirs = bk7.grad_phi(cov)
        require(np.max(np.abs(ours - theirs)) < 2e-12, 'P2 upstream gradient mismatch')
        enclosure = phifun.grad_ad(entries, entries).v
        require(np.all(ours >= enclosure.lo - 2e-13) and np.all(ours <= enclosure.hi + 2e-13),
                'P2 gradient disagrees with upstream point enclosure')
        all_sign_moment = 16 * bk7.Q4(cov) - 1 - sum((2 / math.pi) * math.asin(cov[i, j])
                                                       for i, j in itertools.combinations(range(4), 2))
        moment = fourth_sign_moment(cov, order)
        require(abs(moment - all_sign_moment) < 2e-10, 'fourth sign moment mismatch')
        rows.append({'point': name,
                     'P2_cost_absolute_difference_from_gc4': abs(computed - reference),
                     'P2_gradient_max_difference_from_bk7': float(np.max(np.abs(ours - theirs))),
                     'fourth_moment_difference_from_all_six_Plackett_path': abs(moment - all_sign_moment)})
    block = np.eye(4)
    block[0, 1] = block[1, 0] = 0.6
    block[2, 3] = block[3, 2] = -0.3
    independent = (2 / math.pi) ** 2 * math.asin(0.6) * math.asin(-0.3)
    require(abs(fourth_sign_moment(block, order) - independent) < 2e-14,
            'independent-pair moment identity failed')
    require(abs(fourth_sign_moment(np.eye(4), order)) < 1e-15, 'independent moment not zero')
    return {'upstream_comparisons': rows, 'independent_pair_moment_identity_passed': True}


def check_p4_covariance_derivative(order=64):
    cov = tree_covariance(4)
    grad = cost_gradient(cov, order)
    checks = []
    h = 1e-5
    # Close, distant, and middle boundary pairs probe different Walsh terms.
    for i, j in ((0, 1), (0, 5), (2, 3)):
        direction = np.zeros_like(cov)
        direction[i, j] = direction[j, i] = 1
        minus, plus = cov - h * direction, cov + h * direction
        require(np.linalg.eigvalsh(minus)[0] > 0 and np.linalg.eigvalsh(plus)[0] > 0,
                'finite difference left the positive-definite domain')
        finite = (expected_cost_quadrature(plus, order) - expected_cost_quadrature(minus, order)) / (2 * h)
        error = abs(finite - grad[i, j])
        require(error < 2e-6, 'P4 scalar finite difference disagrees with gradient')
        checks.append({'pair': [i, j], 'analytic_numerical_derivative': float(grad[i, j]),
                       'central_difference': finite, 'absolute_difference': error,
                       'finite_difference_step': h})
    return {'P4_scalar_finite_difference_checks': checks}


def check_cone_chain_rule(order=64):
    k = 4
    types = enumerate_types(k)
    gen = generators(types)
    tree = np.array([float(w) for w in tree_distribution(k, types)])
    # Every coordinate is positive, permitting genuinely two-sided perturbations
    # toward arbitrary generator vertices without leaving the mass simplex.
    masses = 0.9 * tree + 0.1 / len(types)
    masses /= masses.sum()
    cov = covariance(k, masses, gen)
    grad = cost_gradient(cov, order)
    derivatives = cone_directional_derivatives(k, masses, gen, order)
    homogeneity = abs(float(masses @ derivatives))
    require(homogeneity < 2e-12, 'homogeneity identity failed')
    scaled = cone_directional_derivatives(k, 7 * masses, gen, order)
    scale_error = float(np.max(np.abs(scaled - derivatives)))
    require(scale_error < 2e-11, 'mass scale invariance failed')
    support = np.flatnonzero(tree)
    indices = [int(support[0]), int(np.argmin(derivatives)), int(np.argmax(derivatives))]
    h = 1e-6
    rows = []
    for index in dict.fromkeys(indices):
        direction = -masses.copy()
        direction[index] += 1
        plus = covariance(k, masses + h * direction, gen)
        minus = covariance(k, masses - h * direction, gen)
        # Independent finite difference of both covariance normalizations,
        # contracted with the already independently checked correlation gradient.
        finite = float(np.sum(grad * ((plus - minus) / (2 * h))) / k)
        error = abs(finite - derivatives[index])
        require(error < 2e-8, 'double normalization chain rule mismatch')
        rows.append({'type_index': index, 'chain_rule_derivative': float(derivatives[index]),
                     'covariance_finite_difference_derivative': finite, 'absolute_difference': error})
    # A non-symmetric perturbation checks actual covariance and gradient transport
    # under the path reversal, not just equality of identical tree matrices.
    asymmetric = 0.92 * masses
    asymmetric[support[0]] += 0.08
    other_cov = covariance(k, asymmetric, gen)
    permutation = np.arange(k + 2)[::-1]
    other_grad = cost_gradient(other_cov, order)
    permuted_grad = cost_gradient(other_cov[np.ix_(permutation, permutation)], order)
    symmetric_gradient = other_grad + other_grad.T
    transported = np.triu(symmetric_gradient[np.ix_(permutation, permutation)], 1)
    symmetry_error = float(np.max(np.abs(permuted_grad - transported)))
    require(symmetry_error < 2e-11, 'path reversal gradient symmetry failed')
    return {'cone_double_normalization_checks': rows,
            'weighted_directional_derivative_homogeneity_error': homogeneity,
            'mass_scale_invariance_max_error': scale_error,
            'path_reversal_gradient_max_error': symmetry_error}


def check_saved_root(path, order=128):
    record = json.loads(path.read_text())
    require(record['path_vertices'] == 4, 'saved root is not P4')
    require(record['optimizer_success'], 'saved root optimizer did not succeed')
    types = enumerate_types(4)
    gen = generators(types)
    masses = np.zeros(len(types))
    for item in record['support']:
        index = item['type_index']
        require(tuple(item['type']) == types[index], 'saved root type does not match enumeration')
        require(item['candidate_mass'] > 0, 'saved root has nonpositive support mass')
        masses[index] = item['candidate_mass']
    require(abs(float(masses.sum()) - 1) < 1e-12, 'saved root masses do not sum to one')
    deriv = cone_directional_derivatives(4, masses, gen, order)
    residual = float(np.max(np.abs(deriv[masses > 0])))
    require(residual < 1e-8, 'saved root has a large independently recomputed residual')
    return {'saved_candidate_file': path.name, 'recomputed_support_derivative_residual': residual,
            'root_check_quadrature_order': order,
            'interpretation': 'Stationarity diagnostic only; does not establish local or global maximality.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order', type=int, default=64)
    parser.add_argument('--candidate', type=Path,
                        default=Path(__file__).with_name('cone_search_P4.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    require(args.order >= 32, 'verification requires quadrature order at least 32')
    result = {'schema_version': 1,
              'status': 'EXACT FINITE WALSH CHECKS AND NUMERICAL DERIVATIVE DIAGNOSTICS; NOT AN INTERVAL CERTIFICATE',
              'quadrature_order': args.order}
    result.update(check_walsh())
    result.update(check_p2_and_fourth_moments(args.order))
    result.update(check_p4_covariance_derivative(args.order))
    result.update(check_cone_chain_rule(args.order))
    result.update(check_saved_root(args.candidate, 2 * args.order))
    result['all_checks_passed'] = True
    text = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
