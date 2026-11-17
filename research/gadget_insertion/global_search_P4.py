#!/usr/bin/env python3
"""Adversarial multi-start search of the full P4 distance-type relaxation.

Active-set optimization checks numerical derivatives toward EVERY generator.
This is a search for counterexamples and stationary candidates, not an upper
bound certificate. The domain boundary is approached through a positive tree
regularizer; singular points are not certified by this numerical routine.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from path_cone import boundary_owners, enumerate_types, generators, tree_distribution
from path_gradient import cost_gradient, expected_cost_quadrature


class Objective:
    def __init__(self, gen, tree, order=24, regularizer=1e-6):
        self.gen = np.asarray(gen, dtype=float)
        self.k = (self.gen.shape[1]-2)//2
        self.order = order
        self.regularizer = regularizer
        self.tree_gram = (self.gen.T*tree) @ self.gen
        self.incidence = np.zeros((self.k+2, 2*self.k+2))
        for j, owner in enumerate(boundary_owners(self.k)):
            self.incidence[j, owner] = self.incidence[j, self.k+j] = 1

    def covariance_and_adjoint(self, gram):
        diagonal = np.diag(gram)
        roots = np.sqrt(np.outer(diagonal, diagonal))
        vertex = gram/roots
        edge = self.incidence @ vertex @ self.incidence.T
        ed = np.diag(edge)
        er = np.sqrt(np.outer(ed, ed))
        cov = edge/er
        upper = cost_gradient(cov, self.order)/self.k
        adjoint_cov = (upper+upper.T)/2
        adjoint_edge = adjoint_cov/er
        adjoint_edge[np.diag_indices_from(edge)] -= np.sum(adjoint_cov*cov, axis=1)/ed
        adjoint_vertex = self.incidence.T @ adjoint_edge @ self.incidence
        adjoint_gram = adjoint_vertex/roots
        adjoint_gram[np.diag_indices_from(gram)] -= np.sum(adjoint_vertex*vertex, axis=1)/diagonal
        return cov, adjoint_gram

    def gram(self, support, weights):
        g = self.gen[support]
        return (1-self.regularizer)*((g.T*weights) @ g)+self.regularizer*self.tree_gram

    def value_gradient(self, support, weights):
        gram = self.gram(support, weights)
        cov, adjoint = self.covariance_and_adjoint(gram)
        value = expected_cost_quadrature(cov, self.order)/self.k
        g = self.gen[support]
        gradient = (1-self.regularizer)*np.einsum('ti,ij,tj->t', g, adjoint, g)
        return -value, -gradient

    def all_directions(self, support, weights):
        gram = self.gram(support, weights)
        cov, adjoint = self.covariance_and_adjoint(gram)
        scores = (1-self.regularizer)*np.einsum('ti,ij,tj->t', self.gen, adjoint, self.gen)
        return scores-float(weights @ scores[support]), cov


def search(starts=12, seed=73191, order=24, rounds=8, regularizer=1e-6):
    k = 4
    types = enumerate_types(k)
    gen = generators(types)
    tree = np.array([float(x) for x in tree_distribution(k, types)])
    tree_support = np.flatnonzero(tree)
    rng = np.random.default_rng(seed)
    objective = Objective(gen, tree, order, regularizer)
    trials = []
    for trial in range(starts):
        if trial == 0:
            support = tree_support.copy()
            weights = tree[support].copy()
            initialization = 'tree'
        else:
            # All random starts include remote cone generators; half omit the
            # tree face from their active set until column generation finds it.
            count = 8 if trial % 3 else 32
            chosen = rng.choice(len(types), size=count, replace=False)
            support = np.unique(np.r_[tree_support, chosen]) if trial % 2 else chosen
            weights = rng.dirichlet(np.full(len(support), 0.5 if trial % 3 else 2.0))
            initialization = 'sparse_random' if trial % 3 else 'diffuse_random'
        started = time.monotonic()
        history = []
        initial_cov = objective.covariance_and_adjoint(objective.gram(support, weights))[0]
        initial_value = expected_cost_quadrature(initial_cov, order)/k
        for iteration in range(rounds):
            fit = minimize(lambda w: objective.value_gradient(support, w), weights,
                           jac=True, method='SLSQP', bounds=[(0, 1)]*len(support),
                           constraints=[{'type': 'eq', 'fun': lambda w: np.sum(w)-1,
                                         'jac': lambda w: np.ones_like(w)}],
                           options={'ftol': 2e-11, 'maxiter': 160})
            # Repair only last-bit constraint tolerances before forming a PSD
            # Gram matrix; record optimizer status rather than assuming success.
            weights = np.maximum(fit.x, 0)
            weights /= weights.sum()
            keep = weights > 1e-9
            support, weights = support[keep], weights[keep]
            weights /= weights.sum()
            derivatives, cov = objective.all_directions(support, weights)
            maximum = float(derivatives.max())
            history.append({'round': iteration, 'support_size': len(support),
                            'optimizer_success': bool(fit.success),
                            'optimizer_message': str(fit.message),
                            'function_evaluations': int(fit.nfev),
                            'maximum_all_generator_direction': maximum})
            if maximum <= 2e-6:
                break
            # Add several favorable columns, keeping the old feasible point.
            existing = set(support.tolist())
            new = [int(i) for i in np.argsort(derivatives)[::-1]
                   if derivatives[i] > 1e-6 and i not in existing][:12]
            if not new:
                break
            support = np.r_[support, new]
            weights = np.r_[weights, np.zeros(len(new))]
        cov = objective.covariance_and_adjoint(objective.gram(support, weights))[0]
        values = {str(n): expected_cost_quadrature(cov, n)/k for n in (64, 128)}
        total_mass = regularizer*tree
        total_mass[support] += (1-regularizer)*weights
        actual_support = np.flatnonzero(total_mass > 0)
        row = {'trial': trial, 'initialization': initialization,
               'initial_cost_per_vertex': initial_value,
               'final_values_by_quadrature_order': values,
               'history': history,
               'seconds': time.monotonic()-started,
               'mass_support': [{'type_index': int(i), 'type': list(types[i]),
                                 'mass': float(total_mass[i])} for i in actual_support]}
        trials.append(row)
        print(json.dumps({'progress': trial+1, 'starts': starts,
                          'value_order128': values['128'], 'support': len(actual_support),
                          'max_direction': history[-1]['maximum_all_generator_direction'],
                          'seconds': row['seconds']}), flush=True)
    return {'schema_version': 1,
            'status': 'NUMERICAL ADVERSARIAL SEARCH, NOT A GLOBAL UPPER CERTIFICATE',
            'path_vertices': k, 'distance_types': len(types), 'seed': seed,
            'training_quadrature_order': order,
            'regularizer_mass_on_natural_tree': regularizer,
            'all_generators_checked_per_column_round': len(types),
            'interpretation': 'A failure to find a worse point proves no global inequality. All derivatives and quadratures are floating point. Generator directions target (1-epsilon)*e_i+epsilon*m_tree, with the stated regularizer epsilon.',
            'maximum_found_order128': max(r['final_values_by_quadrature_order']['128'] for r in trials),
            'trials': trials}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--starts', type=int, default=12)
    p.add_argument('--seed', type=int, default=73191)
    p.add_argument('--order', type=int, default=24)
    p.add_argument('--rounds', type=int, default=8)
    p.add_argument('--regularizer', type=float, default=1e-6)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    if args.starts < 1 or args.order < 8 or args.rounds < 1 or not 0 < args.regularizer < 1:
        p.error('invalid search settings')
    result = search(args.starts, args.seed, args.order, args.rounds, args.regularizer)
    rendered = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding='utf-8')
    print(rendered, end='')


if __name__ == '__main__':
    main()
