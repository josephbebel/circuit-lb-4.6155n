#!/usr/bin/env python3
"""Numerical covariance Hessian diagnostic; no concavity certificate.

The finite differences below inspect the curvature of Phi_P4(C)/4, with C's
15 upper-triangle entries as coordinates. They also restrict the Hessian to
the tangent plane of the first derivative. Positive numerical curvature in
that plane warns that a supporting linear cut alone is insufficient.
Neither finite-difference agreement nor eigenvalue signs are interval proofs.
"""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from path_cone import enumerate_types, generators, covariance
from path_gradient import cost_gradient


def diagnose(candidate_file, order=48, step=1e-5):
    candidate = json.loads(Path(candidate_file).read_text())
    k = candidate["path_vertices"]
    if k != 4:
        raise ValueError("this diagnostic is specifically for P4")
    types = enumerate_types(k)
    gen = generators(types)
    masses = np.zeros(len(types))
    for item in candidate["support"]:
        index = item["type_index"]
        if tuple(item["type"]) != types[index]:
            raise ValueError("candidate type/index mismatch")
        masses[index] = item["candidate_mass"]
    cov = covariance(k, masses, gen)
    pairs = list(itertools.combinations(range(k+2), 2))
    gradient_matrix = cost_gradient(cov, order)/k
    gradient = np.array([gradient_matrix[i, j] for i, j in pairs])
    _, _, vh = np.linalg.svd(gradient[None, :], full_matrices=True)
    tangent_basis = vh[1:].T
    records = []
    for h in (step, step/2):
        raw = np.empty((len(pairs), len(pairs)))
        for col, (i, j) in enumerate(pairs):
            delta = np.zeros_like(cov)
            delta[i, j] = delta[j, i] = h
            plus = cost_gradient(cov+delta, order)/k
            minus = cost_gradient(cov-delta, order)/k
            raw[:, col] = [(plus[a, b]-minus[a, b])/(2*h) for a, b in pairs]
        hessian = (raw+raw.T)/2
        eigenvalues, eigenvectors = np.linalg.eigh(hessian)
        tangent = tangent_basis.T @ hessian @ tangent_basis
        records.append({
            "difference_step": h,
            "unsymmetrized_max_asymmetry": float(np.max(np.abs(raw-raw.T))),
            "covariance_hessian_eigenvalues": eigenvalues.tolist(),
            "gradient_tangent_hessian_eigenvalues": np.linalg.eigvalsh(tangent).tolist(),
            "positive_covariance_modes": [
                {"eigenvalue": float(value), "coordinates": eigenvectors[:, idx].tolist(),
                 "gradient_dot_direction": float(gradient @ eigenvectors[:, idx])}
                for idx, value in enumerate(eigenvalues) if value > 0],
        })
    return {
        "schema_version": 1,
        "status": "NUMERICAL CURVATURE DIAGNOSTIC; NOT A CERTIFICATE",
        "path_vertices": k,
        "objective": "Phi_P4(C)/4 in upper-triangle covariance coordinates",
        "coordinate_pairs": pairs,
        "point_source": Path(candidate_file).name,
        "quadrature_order": order,
        "covariance_eigenvalues": np.linalg.eigvalsh(cov).tolist(),
        "gradient": gradient.tolist(),
        "finite_differences": records,
        "limitations": [
            "The candidate masses, Gaussian derivatives and finite differences use floating point.",
            "Numerical positive curvature does not by itself certify a strict saddle theorem.",
            "The covariance tangent space is larger than the realizable cone tangent directions.",
            "The diagnostic neither disproves nor proves a proposed global cone bound.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path,
                        default=Path(__file__).with_name("cone_search_P4.json"))
    parser.add_argument("--order", type=int, default=48)
    parser.add_argument("--step", type=float, default=1e-5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.order < 8 or not 0 < args.step <= 1e-3:
        parser.error("require order>=8 and 0<step<=0.001")
    result = diagnose(args.candidate, args.order, args.step)
    rendered = json.dumps(result, indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
