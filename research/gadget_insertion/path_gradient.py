"""Deterministic numerical Gaussian derivatives for P2, P3, and P4.

Boolean Fourier coefficients are computed exactly as integer numerators over a
power of two. Gaussian integrations use floating Gauss-Legendre quadrature.
Increasing the quadrature order is a numerical diagnostic, NOT a certificate.
"""

from functools import lru_cache
import itertools
import math

import numpy as np

from path_gadgets import cost_table


def correlation(gram):
    diagonal = np.diag(gram)
    if np.any(diagonal <= 0):
        raise ValueError("conditional variances must be positive")
    roots = np.sqrt(diagonal)
    return gram / np.outer(roots, roots)


@lru_cache(maxsize=8)
def quadrature(order):
    if not isinstance(order, int) or order < 8:
        raise ValueError("quadrature order must be an integer >= 8")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return (nodes + 1) / 2, weights / 2


def fourth_sign_moment(covariance, order=64):
    """E[sign(X0) sign(X1) sign(X2) sign(X3)], floating approximation.

    Scale only the four cross correlations between pairs (0,1) and (2,3).
    At scale zero, the two pair products are independent. The correlation
    derivative of the fourth sign moment is
       4 asin(rho_conditional) / (pi^2 sqrt(1-r_ij^2)).
    Integrating these four derivatives gives one smooth one-dimensional
    integral. All interpolated matrices are positive definite when C is.
    """
    cov = np.asarray(covariance, dtype=float)
    if cov.shape != (4, 4):
        raise ValueError("a four-dimensional correlation matrix is required")
    nodes, weights = quadrature(order)
    matrices = np.broadcast_to(cov, (order, 4, 4)).copy()
    cross = ((0, 2), (0, 3), (1, 2), (1, 3))
    for i, j in cross:
        matrices[:, i, j] = matrices[:, j, i] = nodes * cov[i, j]
    total = np.zeros(order)
    for i, j in cross:
        k, ell = [x for x in range(4) if x not in (i, j)]
        r = matrices[:, i, j]
        denominator = 1-r*r
        c, d = matrices[:, k, i], matrices[:, ell, i]
        e, f = matrices[:, k, j], matrices[:, ell, j]
        numerator = matrices[:, k, ell]*denominator - (c*d+e*f-r*(c*f+e*d))
        variance_k = denominator-c*c-e*e+2*r*c*e
        variance_l = denominator-d*d-f*f+2*r*d*f
        if np.any(variance_k <= 0) or np.any(variance_l <= 0):
            raise ValueError("quadrature path reached a singular conditional matrix")
        rho = numerator / np.sqrt(variance_k*variance_l)
        if np.any(np.abs(rho) > 1+1e-10):
            raise ValueError("invalid conditional correlation")
        total += cov[i, j]*np.arcsin(np.clip(rho, -1, 1))/np.sqrt(denominator)
    independent = (2/math.pi)**2 * math.asin(cov[0, 1])*math.asin(cov[2, 3])
    return float(independent + (4/math.pi**2)*(weights @ total))


@lru_cache(maxsize=None)
def mixed_difference_fourier(k, i, j):
    """Exact Walsh numerators of Delta_ij F on the remaining coordinates."""
    b = k+2
    rest = tuple(r for r in range(b) if r not in (i, j))
    table = cost_table(k)
    size = 1 << len(rest)
    values = []
    for mask in range(size):
        base = sum(1 << r for pos, r in enumerate(rest) if (mask >> pos) & 1)
        values.append(int(table[base | (1 << i) | (1 << j)])
                      - int(table[base | (1 << i)])
                      - int(table[base | (1 << j)]) + int(table[base]))
    numerators = tuple(sum(values[s] * (-1)**((s & mask).bit_count())
                           for s in range(size)) for mask in range(size))
    if any(numerators[mask] for mask in range(size) if mask.bit_count() % 2):
        raise ArithmeticError("complement symmetry failed")
    return rest, numerators, size


def cost_gradient(covariance, order=64):
    """Derivative of total expected cost w.r.t. each upper-triangle entry."""
    cov = np.asarray(covariance, dtype=float)
    k = len(cov)-2
    if cov.shape != (k+2, k+2) or k not in (2, 3, 4):
        raise ValueError("this derivative evaluator supports P2, P3, and P4")
    gradient = np.zeros_like(cov)
    for i, j in itertools.combinations(range(k+2), 2):
        rest, numerator, denominator = mixed_difference_fourier(k, i, j)
        conditioned = [i, j]
        small = cov[np.ix_(conditioned, conditioned)]
        cross = cov[np.ix_(rest, conditioned)]
        schur = cov[np.ix_(rest, rest)] - cross @ np.linalg.solve(small, cross.T)
        cond = correlation(schur)
        expectation = float(numerator[0])
        for u, v in itertools.combinations(range(k), 2):
            expectation += numerator[(1 << u) | (1 << v)]*(2/math.pi)*math.asin(
                float(np.clip(cond[u, v], -1, 1)))
        if k == 4:
            expectation += numerator[15]*fourth_sign_moment(cond, order)
        expectation /= denominator
        gradient[i, j] = expectation/(2*math.pi*math.sqrt(1-cov[i, j]**2))
    return gradient


def expected_cost_quadrature(covariance, order=64):
    """Total expected cost via the positive-definite path I+t(C-I).

    This has an inner one-dimensional quadrature for P4 and is deterministic.
    No error bound or directed rounding is implemented here.
    """
    cov = np.asarray(covariance, dtype=float)
    k = len(cov)-2
    if k not in (2, 3, 4):
        raise ValueError("point integration supports P2, P3, and P4")
    identity = np.eye(k+2)
    direction = cov-identity
    nodes, weights = quadrature(order)
    value = float(np.mean(cost_table(k)))
    for t, weight in zip(nodes, weights):
        gradient = cost_gradient(identity+t*direction, order)
        value += float(weight * np.sum(direction*gradient))
    return value


def cone_directional_derivatives(k, masses, gen, order=64):
    """Numerical derivatives in all directions m -> (1-t)m+t*e_i.

    Input masses sum to one. Returned derivatives are per vertex. Calculus
    includes both row and edge-score normalizations, including nonunit Gram
    diagonals; no frozen-normalization approximation is used.
    """
    from path_cone import boundary_owners
    masses = np.asarray(masses, dtype=float)
    if masses.shape != (len(gen),) or np.any(masses < 0) or masses.sum() <= 0:
        raise ValueError("invalid masses")
    masses = masses / masses.sum()
    gram = (gen.T*masses) @ gen
    diag = np.diag(gram)
    root_products = np.sqrt(np.outer(diag, diag))
    vertex = gram/root_products
    incidence = np.zeros((k+2, 2*k+2))
    for j, owner in enumerate(boundary_owners(k)):
        incidence[j, owner] = incidence[j, k+j] = 1
    edge = incidence @ vertex @ incidence.T
    edge_diag = np.diag(edge)
    edge_root_products = np.sqrt(np.outer(edge_diag, edge_diag))
    cov = edge/edge_root_products
    gradient = cost_gradient(cov, order)
    delta_gram = gen[:, :, None]*gen[:, None, :]-gram
    ratios = np.diagonal(delta_gram, axis1=1, axis2=2)/diag
    delta_vertex = delta_gram/root_products - vertex*(ratios[:, :, None]+ratios[:, None, :])/2
    delta_edge = np.einsum('iv,tvw,jw->tij', incidence, delta_vertex, incidence)
    edge_ratios = np.diagonal(delta_edge, axis1=1, axis2=2)/edge_diag
    delta_cov = delta_edge/edge_root_products - cov*(edge_ratios[:, :, None]+edge_ratios[:, None, :])/2
    return np.einsum('ij,tij->t', gradient, delta_cov)/k
