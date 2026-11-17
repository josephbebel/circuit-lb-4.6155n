#!/usr/bin/env python3
"""Rigorous trivariate lower bounds obstructing a relaxed H7 window target.

The witness has exact positive rational masses on the ten P8 tree types.
Only square roots, normalizations, and the Gaussian pair-sign arcsine law
are needed. No high-dimensional integration or sampling enters this check.
The conclusion concerns a relaxed covariance face, not graph realizability.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import platform

import mpmath
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
ARITHMETIC_PATH = ROOT / 'mtree' / 'enclosure' / 'ivarith.py'
ARITHMETIC_SHA256 = 'cb1aec517a74efc071374e5cd4426e9bd05958e8d03559ed1d657d6e85186145'
COUNTS = (128, 640, 64, 160, 80, 40, 20, 2, 5, 1)


def require(test, message):
    if not test:
        raise ArithmeticError(message)


require(hashlib.sha256(ARITHMETIC_PATH.read_bytes()).hexdigest() == ARITHMETIC_SHA256,
        'the pinned original interval-arithmetic source has changed')
spec = importlib.util.spec_from_file_location('_window_ivarith', ARITHMETIC_PATH)
arithmetic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(arithmetic)
IA, asin, PI = arithmetic.IA, arithmetic.asin, arithmetic.PI
dn, up = arithmetic.dn, arithmetic.up


def rational_interval(value):
    value = Fraction(value)
    near = float(value)
    actual = Fraction.from_float(near)
    result = IA(dn(near) if actual > value else near,
                up(near) if actual < value else near)
    require(Fraction.from_float(float(result.lo)) <= value <=
            Fraction.from_float(float(result.hi)), 'rational enclosure failed')
    return result


def endpoints(value, places=15):
    lo = Fraction.from_float(float(value.lo))
    hi = Fraction.from_float(float(value.hi))
    require(lo <= hi, 'reversed interval')
    scale = 10 ** places
    lower = (lo.numerator * scale) // lo.denominator
    upper = -((-hi.numerator * scale) // hi.denominator)

    def decimal(n):
        sign = '-' if n < 0 else ''
        n = abs(n)
        return f'{sign}{n // scale}.{n % scale:0{places}d}'

    return {'exact_binary64_endpoint_fractions': [str(lo), str(hi)],
            'outward_decimal_enclosure': [decimal(lower), decimal(upper)]}


def check_pi_enclosure():
    """Check the inherited binary64 pi bounds against an interval constant."""
    require(mpmath.iv.prec >= 80, 'interval-pi check needs at least 80-bit precision')
    lo = Fraction.from_float(float(PI.lo))
    hi = Fraction.from_float(float(PI.hi))
    # Binary64 endpoints are dyadic, with at most 53 significant bits.
    # These conversions are exact at the inherited 80-bit interval precision.
    low = mpmath.iv.mpf(lo.numerator) / mpmath.iv.mpf(lo.denominator)
    high = mpmath.iv.mpf(hi.numerator) / mpmath.iv.mpf(hi.denominator)
    require(bool(low < mpmath.iv.pi) and bool(mpmath.iv.pi < high),
            'the pinned PI interval failed an independent interval-pi check')
    return endpoints(PI)


def tree_types():
    """Derive the ten tree rays directly from integer tree distances."""
    owners = (0, 0, 1, 2, 3, 4, 5, 6, 7, 7)
    locations = [(v, 0) for v in range(8)] + [(v, 1) for v in owners]
    types = []
    for j, owner in enumerate(owners):
        source = 8 + j
        distances = [0 if v == source else abs(location - owner) + extra + 1
                     for v, (location, extra) in enumerate(locations)]
        types.append(tuple(d - distances[0] for d in distances))
    edges = [(v, v + 1) for v in range(7)]
    edges += [(owner, 8 + j) for j, owner in enumerate(owners)]
    for t in types:
        require(t[0] == 0, 'relative-distance normalization failed')
        require(all(abs(t[u] - t[v]) <= 1 for u, v in edges),
                'a type violates edge distance consistency')
        for v in range(8):
            neighbors = [b if a == v else a for a, b in edges if v in (a, b)]
            require(len(neighbors) == 3 and any(t[u] == t[v] - 1 for u in neighbors),
                    'a tree type violates the generic core descent condition')
    require(len(set(types)) == 10, 'tree types are not distinct')
    return owners, types


def covariance_interval(types, masses, owners):
    sqrt2 = IA(2.0).sqrt()

    def power(exponent):
        return rational_interval(Fraction(2) ** (exponent // 2)) * (
            sqrt2 if exponent % 2 else IA(1.0))

    weights = [rational_interval(w) for w in masses]
    gram = [[sum((w * power(-t[i] - t[j]) for t, w in zip(types, weights)), IA(0.0))
             for j in range(18)] for i in range(18)]
    require(all(float(gram[i][i].lo) > 0 for i in range(18)),
            'Gram diagonal positivity failed')
    vertex = [[IA(1.0) if i == j else gram[i][j] / (gram[i][i] * gram[j][j]).sqrt()
               for j in range(18)] for i in range(18)]
    # Boundary indices 1,...,8 are the selected pendant/matching anchors,
    # with respective core owners 0,...,7. The other two are path ends.
    selected = tuple(range(1, 9))
    require(tuple(owners[j] for j in selected) == tuple(range(8)),
            'anchor projection is not one score per core vertex')
    pairs = [(owners[j], 8 + j) for j in selected]
    edge = [[sum((vertex[x][y] for x in pairs[i] for y in pairs[j]), IA(0.0))
             for j in range(8)] for i in range(8)]
    cov = [[IA(1.0) if i == j else edge[i][j] / (edge[i][i] * edge[j][j]).sqrt()
            for j in range(8)] for i in range(8)]
    require(all(-1 < float(cov[i][j].lo) <= float(cov[i][j].hi) < 1
                for i in range(8) for j in range(i)),
            'off-diagonal correlations are not strictly inside (-1,1)')
    return cov


def positive_definiteness(cov):
    """Interval Cholesky with strictly positive pivots."""
    n = len(cov)
    lower = [[IA(0.0) for _ in range(n)] for _ in range(n)]
    pivots = []
    for i in range(n):
        pivot = cov[i][i] - sum((lower[i][k] * lower[i][k] for k in range(i)), IA(0.0))
        require(float(pivot.lo) > 0, 'positive definiteness was not certified')
        pivots.append(endpoints(pivot))
        lower[i][i] = pivot.sqrt()
        for j in range(i + 1, n):
            lower[j][i] = (cov[j][i] - sum((lower[j][k] * lower[i][k]
                                          for k in range(i)), IA(0.0))) / lower[i][i]
    return pivots


def boolean_check():
    for mask in range(256):
        b = [(mask >> i) & 1 for i in range(8)]
        d = [int(b[i] != b[i + 1]) for i in range(7)]

        def window(changes):
            product = 1
            value = 0
            for j in range(7):
                product *= changes[6 - j]
                value += (-1) ** j * product
            return value

        hf, hr = window(d), window(d[::-1])
        af, ar = d[6] * (1 - d[5]), d[0] * (1 - d[1])
        require(hf in (0, 1) and hr in (0, 1), 'H7 range check failed')
        require(hf >= af and hr >= ar and hf + hr >= af + ar,
                'Boolean trivariate lower bound failed')
    return 256


def run():
    pi_enclosure = check_pi_enclosure()
    owners, types = tree_types()
    require(sum(COUNTS) == 1140 and all(x > 0 for x in COUNTS), 'invalid witness counts')
    masses = [Fraction(x, sum(COUNTS)) for x in COUNTS]
    require(sum(masses) == 1, 'rational masses do not sum to one')
    cov = covariance_interval(types, masses, owners)
    pivots = positive_definiteness(cov)
    checked = boolean_check()
    q_forward = IA(0.25) + (asin(cov[5][6]) - asin(cov[6][7]) - asin(cov[5][7])) / (2 * PI)
    q_reverse = IA(0.25) + (asin(cov[1][2]) - asin(cov[0][1]) - asin(cov[0][2])) / (2 * PI)
    q_symmetric = (q_forward + q_reverse) / 2
    frontier = Fraction(1341, 10000)
    clean_lower = Fraction(17, 125)
    bounds = {}
    for name, value in [('forward', q_forward), ('reverse', q_reverse),
                        ('symmetric_average', q_symmetric)]:
        require(Fraction.from_float(float(value.lo)) > clean_lower > frontier,
                f'{name} did not strictly exceed the claimed lower bound')
        bounds[name] = {
            'trivariate_probability': endpoints(value),
            'strict_lower_bound_for_corresponding_H7_expectation': str(clean_lower),
            'excess_over_frontier': endpoints(value - rational_interval(frontier)),
        }
    return {
        'status': 'RIGOROUS RELAXED-CONE OBSTRUCTION; NOT A GRAPH REALIZABILITY CLAIM',
        'path_core_vertices': 8,
        'witness_counts': list(COUNTS), 'witness_count_sum': sum(COUNTS),
        'support': [{'type': list(t), 'mass': str(w)} for t, w in zip(types, masses)],
        'projection_boundary_indices': list(range(1, 9)),
        'projection_owners': list(range(8)),
        'all_boolean_assignments_checked': checked,
        'positive_definiteness_cholesky_pivots': pivots,
        'projected_covariance': [[endpoints(x) for x in row] for row in cov],
        'lower_bounds': bounds,
        'frontier_comparison_target': str(frontier),
        'certified_clean_strict_lower_bound': str(clean_lower),
        'arithmetic': {
            'source': str(ARITHMETIC_PATH.relative_to(ROOT)),
            'sha256': ARITHMETIC_SHA256,
            'method': 'outward binary64 elementary operations; arcsine endpoints verified by mpmath interval sine; exact rational input and output conversion',
            'pi_enclosure_checked_against_mpmath_interval': pi_enclosure,
        },
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                     'mpmath': mpmath.__version__},
        'scope': 'The covariance belongs to the positive ten-type P8 tree face after the anchor projection. This refutes a pointwise H7 or symmetrized-H7 upper bound of 0.1341 on that relaxed face. No cubic graph realizing the masses, no spatial-average obstruction, and no new circuit lower bound is established.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: result[key] for key in (
        'status', 'witness_counts', 'witness_count_sum', 'lower_bounds',
        'certified_clean_strict_lower_bound', 'all_boolean_assignments_checked')}, indent=2))
