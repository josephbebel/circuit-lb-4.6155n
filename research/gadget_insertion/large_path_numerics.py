"""Exploratory conditional-Gaussian calculus for larger path gadgets.

All matrix operations, CDFs and Gaussian expectations are floating point.
There are no rigorous integration errors or global optimization guarantees.
The fixed Cholesky representation used for objective values varies smoothly
with a positive-definite covariance, avoiding eigenbasis jumps in optimization.
"""
from functools import lru_cache
import itertools
import math

import numpy as np
from scipy.special import ndtr

from path_cone import boundary_owners, correlation
from path_gadgets import conditional_values, cost_table


def incidence(k):
    result = np.zeros((k+2, 2*k+2))
    for j, owner in enumerate(boundary_owners(k)):
        result[j, owner] = result[j, k+j] = 1
    return result


def covariance_from_gram(k, gram):
    vertex = correlation(gram)
    matrix = incidence(k)
    return correlation(matrix @ vertex @ matrix.T)


def smooth_cost_samples(covariance, normals):
    """Conditional TOTAL path cost; normals has boundary_size columns.

    If alpha*=1/(1^T C^-1 1), then C-alpha*11^T is positive semidefinite.
    Taking 0.95*alpha* makes the residual positive definite. A common Gaussian
    with coefficient sqrt(alpha) is integrated out, and a Cholesky factor
    represents the residual. This identity is exact; its evaluation is not.
    """
    cov = np.asarray(covariance, dtype=float)
    normals = np.asarray(normals, dtype=float)
    b = len(cov)
    if cov.shape != (b, b) or normals.ndim != 2 or normals.shape[1] != b:
        raise ValueError("incompatible covariance or standard-normal samples")
    ones = np.ones(b)
    alpha = 0.95 / float(ones @ np.linalg.solve(cov, ones))
    if alpha <= 0:
        raise ValueError("positive common-factor variance required")
    residual = cov - alpha * np.ones_like(cov)
    factor = np.linalg.cholesky((residual+residual.T)/2)
    return conditional_values(np.full(b, math.sqrt(alpha)), normals @ factor.T, b-2)


def smooth_expected_cost(covariance, normals):
    return float(np.mean(smooth_cost_samples(covariance, normals)))


@lru_cache(maxsize=None)
def mixed_difference_table(k, i, j):
    rest = tuple(r for r in range(k+2) if r not in (i, j))
    table = cost_table(k)
    result = np.empty(1 << k, dtype=np.int8)
    for mask in range(1 << k):
        base = sum(1 << r for pos, r in enumerate(rest) if mask >> pos & 1)
        result[mask] = int(table[base | (1 << i) | (1 << j)]) - int(table[base | (1 << i)]) - int(table[base | (1 << j)]) + int(table[base])
    return rest, result


def table_conditional_samples(covariance, table, normals):
    """Integrate any Boolean table along one Gaussian factor.

    The common-factor coefficients may have either sign. Before the first
    crossing, negative coefficients have bit 1; every crossing toggles its bit.
    This handles the negative conditional correlations in Plackett gradients.
    """
    cov = np.asarray(covariance, dtype=float)
    table = np.asarray(table, dtype=float)
    d = len(cov)
    normals = np.asarray(normals, dtype=float)
    if normals.ndim != 2 or normals.shape[1] != d-1 or len(table) != 1 << d:
        raise ValueError("incompatible conditional Gaussian table dimensions")
    eig, vec = np.linalg.eigh(cov)
    if eig[0] <= 0:
        raise ValueError("positive-definite conditional covariance required")
    w = vec[:, -1]*math.sqrt(eig[-1])
    residual = vec[:, :-1]*np.sqrt(eig[:-1])[None, :]
    if np.min(np.abs(w)) < 1e-9:
        # A deterministic orthogonal rotation avoids zero common coefficients
        # at symmetric matrices, without changing the represented covariance.
        direction = np.arange(1, d+1, dtype=float)
        direction /= np.linalg.norm(direction)
        v = -direction
        v[0] += 1
        v /= np.linalg.norm(v)
        rotation = np.eye(d)-2*np.outer(v, v)
        factor = np.linalg.cholesky(cov) @ rotation
        w, residual = factor[:, 0], factor[:, 1:]
    if np.any(np.abs(w) < 1e-14):
        raise ValueError("zero common coefficient after deterministic rotation")
    thresholds = -(normals @ residual.T)/w
    order = np.argsort(thresholds, axis=1, kind='stable')
    thresholds = np.take_along_axis(thresholds, order, axis=1)
    initial = sum(1 << i for i in range(d) if w[i] < 0)
    masks = initial ^ np.bitwise_xor.accumulate(np.left_shift(np.uint64(1), order.astype(np.uint64)), axis=1)
    after = np.asarray(table)[masks]
    before = np.column_stack((np.full(len(normals), table[initial]), after[:, :-1]))
    return table[initial] + np.sum((after-before)*ndtr(-thresholds), axis=1)


def cost_gradient_qmc(covariance, normals):
    """Price/Plackett upper-triangle gradient, with numerical conditional means.

    normals has k-1 columns for a Pk gadget. A single sample set is reused for
    all correlation derivatives. Every entry remains an estimate.
    """
    cov = np.asarray(covariance, dtype=float)
    k = len(cov)-2
    grad = np.zeros_like(cov)
    for i, j in itertools.combinations(range(k+2), 2):
        rest, table = mixed_difference_table(k, i, j)
        cross = cov[np.ix_(rest, [i, j])]
        schur = cov[np.ix_(rest, rest)] - cross @ np.linalg.solve(cov[np.ix_([i, j], [i, j])], cross.T)
        conditional = correlation(schur)
        value = float(np.mean(table_conditional_samples(conditional, table, normals)))
        grad[i, j] = value/(2*math.pi*math.sqrt(1-cov[i, j]**2))
    return grad


def correlation_pullback(gram, gradient):
    """Adjoint derivative of corr(gram), for a symmetric matrix gradient."""
    diag = np.diag(gram)
    roots = np.sqrt(np.outer(diag, diag))
    corr = gram/roots
    result = gradient/roots
    result[np.diag_indices_from(result)] -= np.sum(gradient*corr, axis=1)/diag
    return result


def correlation_to_gram_gradient(k, gram, upper_gradient):
    """Pull a correlation gradient back through both normalizations."""
    matrix = incidence(k)
    vertex = correlation(gram)
    edge = matrix @ vertex @ matrix.T
    # Symmetric full-matrix contraction equals an upper-triangle contraction.
    symmetric = (upper_gradient+upper_gradient.T)/2
    vertex_gradient = matrix.T @ correlation_pullback(edge, symmetric) @ matrix
    return correlation_pullback(gram, vertex_gradient)/k


def generator_directions(types, gram, gram_gradient, batch_size=32768):
    """Directional estimates for all types without a giant derivative tensor.

    Each returned value is g_t^T A g_t - <A,G>, for the per-vertex pullback A.
    """
    result = np.empty(len(types))
    baseline = float(np.sum(gram_gradient*gram))
    for start in range(0, len(types), batch_size):
        gen = np.exp2(-np.asarray(types[start:start+batch_size], dtype=float)/2)
        result[start:start+len(gen)] = np.sum((gen @ gram_gradient)*gen, axis=1)-baseline
    return result
