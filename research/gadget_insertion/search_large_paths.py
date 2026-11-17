#!/usr/bin/env python3
"""Adversarial P5--P7 cone experiments with independent numerical holdouts.

The experiment optimizes the tree face, ranks every cone generator using a
numerical gradient, searches sparse mixtures, and re-optimizes selected points.
None of these finite searches supplies an upper bound over the cone.
"""
import argparse
import hashlib
import heapq
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.special import ndtri, softmax
from scipy.stats import qmc

from path_cone import (boundary_owners, compact_types, count_types, covariance,
                       generators, tree_support)
from large_path_numerics import (correlation_to_gram_gradient, cost_gradient_qmc,
                                 generator_directions, smooth_expected_cost)


def normal_points(dimension, power, seed):
    uniform = qmc.Sobol(dimension, scramble=True, seed=seed).random_base2(power)
    return ndtri(np.clip(uniform, np.finfo(float).eps, 1-np.finfo(float).eps))


def candidate(types, masses, label):
    return {'types': np.asarray(types, dtype=np.int8),
            'masses': np.asarray(masses, dtype=float)/np.sum(masses), 'label': label}


def point_covariance(k, point):
    return covariance(k, point['masses'], generators(point['types']))


def point_value(k, point, normals):
    try:
        value = smooth_expected_cost(point_covariance(k, point), normals)/k
        return value if np.isfinite(value) else -np.inf
    except (ValueError, ArithmeticError, np.linalg.LinAlgError):
        return -np.inf


def optimize_point(k, start, normals, maxiter=80):
    gen = generators(start['types'])
    initial = np.maximum(start['masses'], 1e-8)
    initial /= initial.sum()
    invalid = 0

    def objective(logs):
        nonlocal invalid
        masses = softmax(np.r_[logs, 0.0])
        try:
            value = smooth_expected_cost(covariance(k, masses, gen), normals)/k
            if np.isfinite(value):
                return -value
        except (ValueError, ArithmeticError, np.linalg.LinAlgError):
            pass
        invalid += 1
        return 1.0

    fit = minimize(objective, np.log(initial[:-1]/initial[-1]), method='BFGS',
                   options={'maxiter': maxiter, 'gtol': 3e-7})
    weights = softmax(np.r_[fit.x, 0.0])
    point = candidate(start['types'], weights, start['label']+'_optimized')
    point['optimization'] = {'success': bool(fit.success), 'message': str(fit.message),
                             'iterations': int(fit.nit), 'evaluations': int(fit.nfev),
                             'invalid_trial_points': invalid,
                             'training_value': float(-fit.fun)}
    return point


def symmetric_tree_fit(k, support, masses, normals):
    groups = np.array([min(owner, k-1-owner) for owner in boundary_owners(k)])
    count = int(groups.max())+1
    gen = generators(support)

    def weights(logs):
        w = masses*np.exp(np.r_[0.0, logs][groups])
        return w/w.sum()

    def objective(logs):
        return -smooth_expected_cost(covariance(k, weights(logs), gen), normals)/k

    fit = minimize(objective, np.zeros(count-1), method='Nelder-Mead',
                   options={'maxiter': 250, 'xatol': 2e-5, 'fatol': 2e-10})
    point = candidate(support, weights(fit.x), 'symmetric_tree_face')
    point['optimization'] = {'success': bool(fit.success), 'message': str(fit.message),
                             'iterations': int(fit.nit), 'evaluations': int(fit.nfev),
                             'log_multipliers': [0.0]+fit.x.tolist(),
                             'training_value': float(-fit.fun)}
    return point


def public_point(point):
    return {'label': point['label'], 'support': [
        {'type': list(map(int, row)), 'mass': float(weight)}
        for row, weight in zip(point['types'], point['masses'])],
        **({'optimization': point['optimization']} if 'optimization' in point else {})}


def search(k, train_power=15, holdout_power=19, repeats=6, random_trials=2000,
           seed=37001, polish=4):
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    types = compact_types(k)
    support, exact = tree_support(k)
    masses = np.array([float(w) for w in exact])
    natural = candidate(support, masses, 'natural_tree')
    training = normal_points(k+2, train_power, seed+1)
    screening = training[:min(len(training), 1024)]
    symmetric = symmetric_tree_fit(k, support, masses, training)
    face_starts = [symmetric]
    for scale in (0.08, 0.30):
        perturbed = symmetric['masses']*np.exp(rng.normal(0, scale, len(masses)))
        face_starts.append(candidate(support, perturbed, f'tree_face_random_start_{scale}'))
    face_points = [symmetric]+[optimize_point(k, start, training) for start in face_starts]
    best_face = max(face_points, key=lambda point: point_value(k, point, training))
    best_face['label'] = 'best_tree_face'
    print(json.dumps({'stage': 'tree_face', 'k': k,
                      'natural_training': point_value(k, natural, training),
                      'best_training': point_value(k, best_face, training)}), flush=True)

    cov = point_covariance(k, best_face)
    gradient_normals = normal_points(k-1, 16, seed+2)
    grad = cost_gradient_qmc(cov, gradient_normals)
    gen = generators(best_face['types'])
    gram = (gen.T*best_face['masses']) @ gen
    pullback = correlation_to_gram_gradient(k, gram, grad)
    directions = generator_directions(types, gram, pullback)
    is_tree = np.zeros(len(types), dtype=bool)
    for row in support:
        is_tree |= np.all(types == np.asarray(row), axis=1)
    off_indices = np.flatnonzero(~is_tree)
    ranked = off_indices[np.argsort(directions[off_indices])[-256:][::-1]]
    print(json.dumps({'stage': 'all_generator_screen', 'k': k, 'types': len(types),
                      'largest_off_face_derivative_estimate': float(directions[ranked[0]])}), flush=True)

    heap = []
    serial = 0
    invalid = 0

    def retain(point):
        nonlocal serial, invalid
        score = point_value(k, point, screening)
        serial += 1
        if not np.isfinite(score):
            invalid += 1
            return
        entry = (score, serial, point)
        if len(heap) < 24:
            heapq.heappush(heap, entry)
        elif score > heap[0][0]:
            heapq.heapreplace(heap, entry)

    for index in ranked[:24]:
        for beta in (0.01, 0.03, 0.1, 0.25, 0.5, 0.8):
            retain(candidate(np.vstack((best_face['types'], types[index])),
                             np.r_[(1-beta)*best_face['masses'], beta],
                             f'single_off_face_{int(index)}_mass_{beta}'))
    for trial in range(random_trials):
        count = int(rng.integers(2, 7))
        pool = ranked if trial % 2 else off_indices
        chosen = rng.choice(pool, size=count, replace=False)
        beta = float(rng.uniform(0.02, 0.98))
        extra = rng.dirichlet(np.full(count, 0.5))
        retain(candidate(np.vstack((best_face['types'], types[chosen])),
                         np.r_[(1-beta)*best_face['masses'], beta*extra],
                         f'random_sparse_{trial}'))
    shortlist = [entry[2] for entry in heap]
    shortlist.sort(key=lambda point: point_value(k, point, training), reverse=True)
    polished = [optimize_point(k, point, training, maxiter=100) for point in shortlist[:polish]]
    all_sparse = shortlist[:polish]+polished
    all_sparse.sort(key=lambda point: point_value(k, point, training), reverse=True)
    validated = [natural, best_face]+all_sparse[:3]
    # Retain separate initializations in the holdout even when optimization
    # brings them very close. This reveals rather than hides such convergence.
    print(json.dumps({'stage': 'sparse_search', 'k': k, 'tested': serial,
                      'best_sparse_training': point_value(k, all_sparse[0], training)}), flush=True)

    covariances = [point_covariance(k, point) for point in validated]
    holdouts = []
    batch_power = min(15, holdout_power)
    for repeat in range(repeats):
        sampler = qmc.Sobol(k+2, scramble=True, seed=seed+1000+repeat)
        totals = np.zeros(len(validated))
        for _ in range(2**(holdout_power-batch_power)):
            uniform = sampler.random(2**batch_power)
            normals = ndtri(np.clip(uniform, np.finfo(float).eps, 1-np.finfo(float).eps))
            for j, point_cov in enumerate(covariances):
                totals[j] += smooth_expected_cost(point_cov, normals)*len(normals)/k
        holdouts.append((totals/2**holdout_power).tolist())
    values = np.array(holdouts)
    records = []
    for j, point in enumerate(validated):
        mean = float(np.mean(values[:, j]))
        difference = values[:, j]-values[:, 0]
        records.append({**public_point(point),
                        'training_cost_per_vertex': point_value(k, point, training),
                        'holdout_mean_cost_per_vertex': mean,
                        'holdout_replicates': values[:, j].tolist(),
                        'holdout_standard_error_heuristic': float(np.std(values[:, j], ddof=1)/np.sqrt(repeats)),
                        'paired_gain_over_tree_mean': float(difference.mean()),
                        'paired_gain_standard_error_heuristic': float(np.std(difference, ddof=1)/np.sqrt(repeats)),
                        'target_minus_mean': {'0.13614': 0.13614-mean, '0.13615': 0.13615-mean}})
    return {
        'schema_version': 1,
        'status': 'ADVERSARIAL NUMERICAL SEARCH; NO RIGOROUS UPPER BOUND OR CIRCUIT RECORD',
        'path_vertices': k, 'counts': count_types(k),
        'sorted_int8_type_sha256': hashlib.sha256(types.tobytes()).hexdigest(),
        'compact_type_bytes': types.nbytes,
        'training': {'seed': seed+1, 'points': 2**train_power,
                     'smooth_common_factor_fraction': 0.95,
                     'tree_face_initializations': len(face_starts)+1},
        'generator_direction_screen': {
            'seed': seed+2, 'points': 2**16, 'all_generators_evaluated': len(types),
            'conditional_gradient_status': 'QMC estimate, no error guarantee',
            'largest_off_face_directions': [
                {'type_index': int(i), 'type': list(map(int, types[i])), 'derivative_estimate': float(directions[i])}
                for i in ranked[:12]]},
        'sparse_search': {'single_type_grid_trials': 24*6, 'random_sparse_trials': random_trials,
                          'coarse_points': len(screening), 'invalid_trial_points': invalid,
                          'shortlist_size': len(shortlist), 'polished_points': len(polished)},
        'holdout': {'points_per_scramble': 2**holdout_power, 'scrambles': repeats,
                    'seeds': list(range(seed+1000, seed+1000+repeats)),
                    'uncertainty': 'Scramble standard errors are heuristic, not rigorous integration errors or confidence guarantees.'},
        'validated_candidates': records,
        'interpretation': [
            'Failure to find a point above a target is not a proof that none exists.',
            'The all-generator scan is a numerical first-derivative calculation at one point, not a global bound.',
            'The finite random and local searches do not exhaust sparse mixtures or the cone.',
            'A cone point need not correspond to a realized cubic graph.',
            'Every claimed universal improvement still needs a full rigorous Gaussian inequality.'
        ],
        'wall_seconds': time.perf_counter()-started,
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, choices=[5, 6, 7], required=True)
    parser.add_argument('--train-power', type=int, default=15)
    parser.add_argument('--holdout-power', type=int, default=19)
    parser.add_argument('--repeats', type=int, default=6)
    parser.add_argument('--random-trials', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=37001)
    parser.add_argument('--polish', type=int, default=4)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 10 <= args.train_power <= 19 or not 12 <= args.holdout_power <= 22:
        parser.error('training power must be10..19 and holdout power12..22')
    if not 2 <= args.repeats <= 32 or not 0 <= args.random_trials <= 50000 or not 1 <= args.polish <= 12:
        parser.error('invalid repeats, random trials, or polish count')
    result = search(args.size, args.train_power, args.holdout_power, args.repeats,
                    args.random_trials, args.seed, args.polish)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'stage': 'complete', 'path': str(args.output), 'k': args.size,
                      'best_holdout_mean': max(x['holdout_mean_cost_per_vertex'] for x in result['validated_candidates']),
                      'wall_seconds': result['wall_seconds']}), flush=True)


if __name__ == '__main__':
    main()
