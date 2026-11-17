"""Exact Boolean path gadgets and exploratory Gaussian integration.

None of the numerical integration routines is an interval certificate.  In
particular, a tree value is not an upper bound over all cubic neighborhoods.
Only NumPy and SciPy are required beyond the Python standard library.

Coordinates: a k-vertex path has internal edges (j,j+1), j=0,...,k-2.
Its k+2 boundary incidences have owners [0,0,1,...,k-2,k-1,k-1].
Boolean 1 means an edge lies below the current threshold. Complementing all
bits leaves the cost unchanged.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math

import numpy as np
from scipy.special import ndtr


def boundary_owners(path_size):
    if not isinstance(path_size, int) or path_size < 2:
        raise ValueError("path_size must be an integer at least 2")
    return [0, 0] + list(range(1, path_size - 1)) + [path_size - 1] * 2


def _sign_array(signs):
    bits = np.asarray(signs)
    if bits.ndim != 2 or bits.shape[1] < 4:
        raise ValueError("signs must have shape (samples, path_size + 2)")
    if not np.all((bits == 0) | (bits == 1)):
        raise ValueError("signs must contain only Boolean values or 0/1")
    return bits.astype(bool, copy=False)


def cost(signs):
    """Return the minimum split-vertex count, by an exact integer O(k) DP."""
    bits = _sign_array(signs)
    n = bits.shape[1] - 2
    d0 = (bits[:, 0] | bits[:, 1]).astype(np.int32)
    d1 = (~bits[:, 0] | ~bits[:, 1]).astype(np.int32)
    for v in range(1, n - 1):
        b = bits[:, v + 1]
        z = np.minimum(d0 + b, d1 + 1)
        o = np.minimum(d0 + 1, d1 + ~b)
        d0, d1 = z, o
    c0 = d0 + (bits[:, -2] | bits[:, -1])
    c1 = d1 + (~bits[:, -2] | ~bits[:, -1])
    return np.minimum(c0, c1)


def least_minimizer(signs):
    """Return the coordinatewise least optimal internal assignment.

    The forward DP stores predecessor choices, using zero to break ties.
    Backtracking selects the lexicographically least reversed bit word. Since
    the minimizers of this submodular function are closed under intersection,
    their coordinatewise least element exists and is least in this order too.
    The separate exhaustive checker verifies this property for small paths.
    """
    bits = _sign_array(signs)
    n = bits.shape[1] - 2
    samples = len(bits)
    d0 = (bits[:, 0] | bits[:, 1]).astype(np.int32)
    d1 = (~bits[:, 0] | ~bits[:, 1]).astype(np.int32)
    predecessors = []
    for v in range(1, n - 1):
        b = bits[:, v + 1]
        p0 = d1 + 1 < d0 + b
        p1 = d1 + ~b < d0 + 1
        predecessors.append(np.column_stack((p0, p1)))
        d0, d1 = np.minimum(d0 + b, d1 + 1), np.minimum(d0 + 1, d1 + ~b)
    out = np.empty((samples, n - 1), dtype=bool)
    out[:, -1] = d1 + (~bits[:, -2] | ~bits[:, -1]) < d0 + (bits[:, -2] | bits[:, -1])
    row = np.arange(samples)
    for j in range(n - 3, -1, -1):
        out[:, j] = predecessors[j][row, out[:, j + 1].astype(int)]
    return out


@lru_cache(maxsize=None)
def cost_table(path_size):
    """Exact table indexed by little-endian boundary bit mask (k <= 18)."""
    boundary_owners(path_size)
    if path_size > 18:
        raise ValueError("the exhaustive cost table is limited to path_size <= 18")
    b = path_size + 2
    masks = np.arange(1 << b, dtype=np.uint64)
    bits = (masks[:, None] >> np.arange(b, dtype=np.uint64)) & 1
    result = cost(bits)
    result.flags.writeable = False
    return result


@dataclass(frozen=True)
class Qsqrt2:
    """An exact a+b*sqrt(2), with rational coefficients."""
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @staticmethod
    def coerce(other):
        return other if isinstance(other, Qsqrt2) else Qsqrt2(Fraction(other))

    def __add__(self, other):
        other = self.coerce(other)
        return Qsqrt2(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __mul__(self, other):
        other = self.coerce(other)
        return Qsqrt2(self.a * other.a + 2 * self.b * other.b,
                      self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        denominator = other.a ** 2 - 2 * other.b ** 2
        if denominator == 0:
            raise ZeroDivisionError
        return self * Qsqrt2(other.a / denominator, -other.b / denominator)

    def __float__(self):
        return float(self.a) + float(self.b) * math.sqrt(2)

    def json(self):
        return {"rational": str(self.a), "sqrt2_coefficient": str(self.b)}


def exact_vertex_corr(distance):
    """r_d=(1+d/3)*2^(-d/2), in Q(sqrt(2))."""
    if not isinstance(distance, int) or distance < 0:
        raise ValueError("distance must be a nonnegative integer")
    coefficient = Fraction(distance + 3, 3 * 2 ** ((distance + 1) // 2))
    return Qsqrt2(coefficient, Fraction(0)) if distance % 2 == 0 else Qsqrt2(Fraction(0), coefficient)


def vertex_corr(distance):
    return (1 + distance / 3) * 2.0 ** (-distance / 2)


def exact_tree_covariance(path_size):
    owners = boundary_owners(path_size)
    denominator = 2 + 2 * exact_vertex_corr(1)
    result = []
    for i, u in enumerate(owners):
        row = []
        for j, v in enumerate(owners):
            d = abs(u - v)
            row.append(Qsqrt2(Fraction(1)) if i == j else
                       (exact_vertex_corr(d) + 2 * exact_vertex_corr(d + 1)
                        + exact_vertex_corr(d + 2)) / denominator)
        result.append(row)
    return result


def tree_covariance(path_size):
    return np.array([[float(x) for x in row] for row in exact_tree_covariance(path_size)])


def gaussian_factor(covariance, tolerance=1e-10):
    """Perron common factor plus residual factors; floating-point only.

    Columns correspond to descending eigenvalues. Small eigenvalues caused
    by roundoff are clipped to zero. This is not a PSD certificate.
    """
    cov = np.asarray(covariance, dtype=float)
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1] or len(cov) < 2:
        raise ValueError("covariance must be square")
    if not np.isfinite(cov).all() or not np.allclose(cov, cov.T, atol=1e-12, rtol=0):
        raise ValueError("covariance must be finite and symmetric")
    if not np.allclose(np.diag(cov), 1, atol=1e-10, rtol=0):
        raise ValueError("covariance must be a correlation matrix")
    eig, vec = np.linalg.eigh(cov)
    if eig[0] < -tolerance:
        raise ValueError("covariance has a negative eigenvalue")
    eig, vec = eig[::-1], vec[:, ::-1]
    leading = vec[:, 0] * math.sqrt(max(eig[0], 0))
    if leading.sum() < 0:
        leading = -leading
    if np.any(leading <= 0):
        raise ValueError("a strictly positive common Gaussian factor is required")
    residual = vec[:, 1:] * np.sqrt(np.maximum(eig[1:], 0))[None, :]
    return leading, residual, eig


def conditional_values(leading, residual_values, path_size):
    """Integrate the remaining one-dimensional Gaussian exactly in formula.

    Arithmetic and normal CDF evaluations use ordinary floating point.
    Returned values are conditional expected TOTAL costs per residual sample.
    """
    leading = np.asarray(leading, dtype=float)
    residual_values = np.asarray(residual_values, dtype=float)
    b = path_size + 2
    if leading.shape != (b,) or residual_values.ndim != 2 or residual_values.shape[1] != b:
        raise ValueError("incompatible leading factor or residual array")
    if not np.all(leading > 0):
        raise ValueError("the common factor must be positive in every coordinate")
    thresholds = -residual_values / leading
    order = np.argsort(thresholds, axis=1, kind="stable")
    sorted_thresholds = np.take_along_axis(thresholds, order, axis=1)
    masks = np.bitwise_or.accumulate(np.left_shift(np.uint64(1), order.astype(np.uint64)), axis=1)
    interval_costs = cost_table(path_size)[masks[:, :-1]]
    interval_masses = np.diff(ndtr(sorted_thresholds), axis=1)
    # Before the first and after the last threshold every bit agrees: cost zero.
    return np.sum(interval_costs * interval_masses, axis=1)


class ConditionalGaussianCost:
    """Cached numerical evaluator of a path cost for an arbitrary covariance."""
    def __init__(self, covariance, path_size=None):
        covariance = np.asarray(covariance, dtype=float)
        self.path_size = len(covariance) - 2 if path_size is None else path_size
        if len(covariance) != self.path_size + 2:
            raise ValueError("covariance dimension does not match path_size")
        self.leading, self.residual, self.eigenvalues = gaussian_factor(covariance)
        cost_table(self.path_size)

    def sample_values(self, normals):
        normals = np.asarray(normals, dtype=float)
        if normals.ndim != 2 or normals.shape[1] != self.residual.shape[1]:
            raise ValueError("normals must have shape (samples, boundary_size - 1)")
        return conditional_values(self.leading, normals @ self.residual.T, self.path_size)


def expected_cost(covariance, normals, path_size=None):
    """Heuristic numerical expected TOTAL cost; never an interval certificate."""
    return float(np.mean(ConditionalGaussianCost(covariance, path_size).sample_values(normals)))


def common_tree_covariance(max_path_size):
    """One tree process containing all prefix paths up to max_path_size.

    Edges: left exterior edge; for each backbone vertex i, its pendant edge
    and the forward backbone edge i--i+1. The last forward edge is exterior.
    A vertex is represented by (backbone_position, is_pendant).
    """
    boundary_owners(max_path_size)
    edges = [((-1, False), (0, False))]
    for i in range(max_path_size):
        edges.extend([((i, False), (i, True)), ((i, False), (i + 1, False))])

    def dist(u, v):
        return 0 if u == v else abs(u[0] - v[0]) + int(u[1]) + int(v[1])

    denom = 2 + 2 * vertex_corr(1)
    cov = np.empty((len(edges), len(edges)))
    for i, edge in enumerate(edges):
        for j, other in enumerate(edges):
            cov[i, j] = sum(vertex_corr(dist(u, v)) for u in edge for v in other) / denom
    np.fill_diagonal(cov, 1.0)
    return cov


def prefix_indices(path_size):
    boundary_owners(path_size)
    return [0] + [1 + 2 * i for i in range(path_size)] + [2 * path_size]
