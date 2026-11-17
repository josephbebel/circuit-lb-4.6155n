#!/usr/bin/env python3
"""Reproducible conditional QMC study; no universal bound or certification.

All requested paths are prefix paths in one common Gaussian tree process.
Each independent Sobol scramble supplies common random numbers to all paths.
The global positive leading Gaussian factor is integrated analytically by
sorting its crossing thresholds; the residual Gaussian integral is estimated.
"""
import argparse
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.special import ndtri
from scipy.stats import qmc

from path_gadgets import (common_tree_covariance, conditional_values,
                          exact_tree_covariance, gaussian_factor, prefix_indices)


TREE_BASELINE = 0.276583606882373471022763934878 / 2


def study(sizes, power=18, repeats=4, seed=1729, batch_power=14):
    sizes = sorted(set(sizes))
    common = common_tree_covariance(max(sizes))
    leading, residual, eig = gaussian_factor(common)
    # Exact zero modes arise from tree-wave identities. Removing modes below
    # a floating tolerance is a numerical efficiency choice, not certification.
    selected = eig[1:] > 1e-12
    residual = residual[:, selected]
    paths = {n: prefix_indices(n) for n in sizes}
    estimates = {n: [] for n in sizes}
    use_batch = min(batch_power, power)
    for repeat in range(repeats):
        sampler = qmc.Sobol(residual.shape[1], scramble=True, seed=seed + repeat)
        totals = {n: 0.0 for n in sizes}
        for _ in range(2 ** (power - use_batch)):
            uniform = sampler.random(2 ** use_batch)
            # Scrambled Sobol points are interior in normal operation. Clipping
            # prevents infinities if a library ever supplies an endpoint.
            uniform = np.clip(uniform, np.finfo(float).eps, 1 - np.finfo(float).eps)
            residual_samples = ndtri(uniform) @ residual.T
            for n, idx in paths.items():
                values = conditional_values(leading[idx], residual_samples[:, idx], n)
                totals[n] += float(np.sum(values))
        for n in sizes:
            estimates[n].append(totals[n] / (2 ** power * n))
    rows = []
    for n in sizes:
        vals = np.asarray(estimates[n])
        mean = float(np.mean(vals))
        row = {"path_vertices": n, "boundary_edges": n + 2,
               "mean_tree_cost_per_vertex": mean,
               "replicate_values": vals.tolist(),
               "replicate_min": float(vals.min()), "replicate_max": float(vals.max()),
               "standard_error_across_scrambles_heuristic": float(np.std(vals, ddof=1) / np.sqrt(repeats)),
               "tree_value_transformed_as_if_universal": 1 + 1 / (2 * mean),
               "transformed_value_status": "CONDITIONAL ILLUSTRATION ONLY; neither a proved bound nor a record"}
        if 2 in estimates:
            diff = np.asarray(estimates[2]) - vals
            row["paired_improvement_over_P2"] = float(diff.mean())
            row["paired_improvement_standard_error_heuristic"] = float(np.std(diff, ddof=1) / np.sqrt(repeats))
        if n == 2:
            row["known_tree_reference_numerical"] = TREE_BASELINE
            row["difference_from_reference"] = mean - TREE_BASELINE
        rows.append(row)
    return {"schema_version": 1,
            "status": "EXPLORATORY NUMERICS, NOT A CERTIFICATE OR UNIVERSAL CUTWIDTH BOUND",
            "method": "common-tree conditional scrambled Sobol; exact Boolean cost lookup; floating Gaussian CDF",
            "points_per_scramble": 2 ** power, "scrambles": repeats,
            "first_seed": seed, "seeds": list(range(seed, seed + repeats)),
            "batch_power": use_batch, "common_covariance_dimension": len(common),
            "residual_gaussian_dimension": residual.shape[1],
            "uncertainty": "Variation across scrambles is heuristic and is not a rigorous error bound or confidence guarantee.",
            "scope": "Tree marginals only, not worst-neighborhood bounds. The L=4 triangle-free kernel scheme requires uniform certification of P4 through P7; its deterministic packing argument is in PROOF.md.",
            "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
            "upstream_commit": "a56500014ba159397903dea81d41ac836f6880d5",
            "results": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes', type=int, nargs='+', default=[2, 3, 4, 6, 8, 12])
    parser.add_argument('--power', type=int, default=18)
    parser.add_argument('--repeats', type=int, default=4)
    parser.add_argument('--seed', type=int, default=1729)
    parser.add_argument('--batch-power', type=int, default=14)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--exact-covariances', type=Path,
                        help='also write exact coefficients a+b*sqrt(2) for each tree covariance')
    args = parser.parse_args()
    if not args.sizes or any(n < 2 or n > 18 for n in args.sizes):
        parser.error('path sizes must be in [2, 18]')
    if not 4 <= args.power <= 26 or not 2 <= args.repeats <= 64:
        parser.error('--power must be in [4,26] and --repeats in [2,64]')
    if not 4 <= args.batch_power <= 20:
        parser.error('--batch-power must be in [4,20]')
    result = study(args.sizes, args.power, args.repeats, args.seed, args.batch_power)
    rendered = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')
    if args.exact_covariances:
        exact = {str(n): [[x.json() for x in row] for row in exact_tree_covariance(n)] for n in sorted(set(args.sizes))}
        args.exact_covariances.parent.mkdir(parents=True, exist_ok=True)
        args.exact_covariances.write_text(json.dumps({"representation": "rational + sqrt2_coefficient * sqrt(2)", "matrices": exact}, indent=2) + '\n')


if __name__ == '__main__':
    main()
