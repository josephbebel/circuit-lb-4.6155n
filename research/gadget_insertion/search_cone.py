#!/usr/bin/env python3
"""Locate and diagnose a tree-face stationary point; exploration, not proof.

The optimizer uses deterministic numerical Plackett derivatives. The point
values are evaluated at two quadrature orders and with independent conditional
Sobol samples. No sampled maximum, gradient test, or numerical agreement is a
universal cone certificate or a rigorous enclosure.
"""

import argparse
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.optimize import root
from scipy.special import ndtri, softmax
from scipy.stats import qmc

from path_cone import enumerate_types, generators, tree_distribution, covariance
from path_gradient import cone_directional_derivatives, expected_cost_quadrature
from path_gadgets import expected_cost


def search(k=4, order=64, power=18, repeats=4, seed=91021):
    types = enumerate_types(k)
    gen = generators(types)
    tree = np.array([float(w) for w in tree_distribution(k, types)])
    support = np.flatnonzero(tree)
    tree_face = gen[support]
    initial = tree[support]

    def weights(log_ratios):
        return softmax(np.r_[log_ratios, 0.0])

    def stationarity(log_ratios):
        w = weights(log_ratios)
        derivative = cone_directional_derivatives(k, w, tree_face, order)
        return derivative[:-1]-derivative[-1]

    fit = root(stationarity, np.log(initial[:-1]/initial[-1]), method='hybr',
               options={'xtol': 1e-9})
    candidate_face = weights(fit.x)
    candidate = np.zeros_like(tree)
    candidate[support] = candidate_face
    cov_tree = covariance(k, tree, gen)
    cov_candidate = covariance(k, candidate, gen)
    numerical_values = {}
    for n in (order, 2*order):
        numerical_values[str(n)] = {
            "tree_cost_per_vertex": expected_cost_quadrature(cov_tree, n)/k,
            "candidate_cost_per_vertex": expected_cost_quadrature(cov_candidate, n)/k,
        }
    independent = []
    for repeat in range(repeats):
        uniforms = qmc.Sobol(k+1, scramble=True, seed=seed+repeat).random_base2(power)
        uniforms = np.clip(uniforms, np.finfo(float).eps, 1-np.finfo(float).eps)
        normals = ndtri(uniforms)
        a = expected_cost(cov_tree, normals)/k
        b = expected_cost(cov_candidate, normals)/k
        independent.append({"seed": seed+repeat, "tree": a, "candidate": b,
                            "paired_candidate_minus_tree": b-a})
    at_tree = cone_directional_derivatives(k, tree, gen, 2*order)
    at_candidate = cone_directional_derivatives(k, candidate, gen, 2*order)
    off_face = tree == 0
    top = np.argsort(at_candidate)[-10:][::-1]
    return {
        "schema_version": 1,
        "status": "NUMERICAL RESEARCH DIAGNOSTIC; NO RIGOROUS UPPER BOUND OR NEW RECORD",
        "path_vertices": k,
        "distance_types": len(types),
        "optimizer": "root of tree-face numerical stationarity equations in log mass ratios",
        "optimizer_success": bool(fit.success),
        "optimizer_message": str(fit.message),
        "stationarity_max_residual": float(np.max(np.abs(stationarity(fit.x)))),
        "support": [{"type_index": int(i), "type": list(types[i]),
                     "tree_mass": float(tree[i]), "candidate_mass": float(candidate[i])}
                    for i in support],
        "deterministic_quadrature_orders": numerical_values,
        "conditional_sobol_validation": {
            "points_per_scramble": 2**power, "scrambles": repeats,
            "values": independent,
            "uncertainty": "Paired estimates and quadrature-order agreement are heuristic, not certified error bounds."
        },
        "directional_diagnostics": {
            "direction": "m -> (1-t)m + t*e_i; derivative of expected cost per vertex at t=0",
            "tree_max_off_face_direction": float(np.max(at_tree[off_face])),
            "tree_max_all_directions": float(np.max(at_tree)),
            "candidate_max_off_face_direction": float(np.max(at_candidate[off_face])),
            "candidate_max_support_absolute_derivative": float(np.max(np.abs(at_candidate[support]))),
            "candidate_largest_directions": [{"type_index": int(i), "type": list(types[i]),
                                              "derivative": float(at_candidate[i])} for i in top],
        },
        "interpretation": [
            "A stationary point need not be a local or global maximum.",
            "Negative outward first derivatives do not establish global optimality.",
            "A relaxed cone point need not be realized by an actual graph.",
            "A value above the tree point would obstruct a tree-extremality proof for this relaxation, but needs a rigorous comparison to become a theorem.",
            "No graph-wide circuit or cutwidth bound is established here."
        ],
        "upstream_commit": "a56500014ba159397903dea81d41ac836f6880d5",
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, choices=[2, 3, 4], default=4)
    parser.add_argument('--order', type=int, default=64)
    parser.add_argument('--power', type=int, default=18)
    parser.add_argument('--repeats', type=int, default=4)
    parser.add_argument('--seed', type=int, default=91021)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.order < 8 or not 4 <= args.power <= 24 or not 1 <= args.repeats <= 32:
        parser.error('require order >= 8, power in [4,24], repeats in [1,32]')
    result = search(args.size, args.order, args.power, args.repeats, args.seed)
    rendered = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding='utf-8')
    print(rendered, end='')


if __name__ == '__main__':
    main()
