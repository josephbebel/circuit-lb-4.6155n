#!/usr/bin/env python3
"""Certify a wide one-sided natural-P4 neighborhood by monotonicity.

The only Gaussian integration is at C_tree-r*(J-I). Every correlation
matrix above it entrywise has no larger submodular sign cost. This does
not assert a global bound on the distance-type cone.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

from p4_interval import (IA, Jet, PI, INDEPENDENT_COST, coefficients,
                         exact_endpoints, integrate, interval_fraction,
                         load_candidate, prove_pd, require)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius', default='1/300')
    parser.add_argument('--tolerance', type=float, default=0.000004)
    parser.add_argument('--degree', type=int, default=4)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    require(args.degree in (4,6,8), 'unsupported Taylor degree')
    radius = F(args.radius)
    require(0 < radius <= F(1,100), 'comparison radius outside checked range')
    directory = Path(__file__).resolve().parent
    repository = directory.parents[2]
    arithmetic_path = repository/'mtree'/'enclosure'/'ivarith.py'
    arithmetic_bytes = arithmetic_path.read_bytes()
    arithmetic_blob = hashlib.sha1(b'blob '+str(len(arithmetic_bytes)).encode()+b'\0'+arithmetic_bytes).hexdigest()
    require(arithmetic_blob=='b99d01f757f2b3d75151ba9403a547c6b8c60477',
            'interval arithmetic differs from the audited upstream version')
    from mpmath import iv
    require(bool(iv.mpf(float(PI.lo))<iv.pi) and bool(iv.pi<iv.mpf(float(PI.hi))),
            'pi interval failed its independent high-precision check')
    sys.path.insert(0,str(directory.parent))
    from path_cone import enumerate_types
    types = enumerate_types(4)
    require(len(types)==4655 and all(t[0]==0 and max(map(abs,t))<=4 for t in types),
            'P4 generator-range check failed')
    candidate_path = directory.parent/'cone_search_P4.json'
    cov, masses = load_candidate(candidate_path, natural_tree=True)
    require(all(sum((F(row['mass'])*F(2)**(-row['type'][i]) for row in masses),F(0))==1
                for i in range(10)), 'natural Gram diagonals are not all exactly one')
    radius_iv = interval_fraction(radius)
    comparison = [[IA(1.0) if i==j else cov[i][j]-radius_iv
                   for j in range(6)] for i in range(6)]
    pivots = prove_pd(comparison)
    full_box_pivots = prove_pd(cov,interval_fraction(5*radius))
    Jet.configure(args.degree)
    integral, stats = integrate(comparison,args.tolerance,initial=8,
                                batch_size=64,max_cells=100000)
    total_cost = interval_fraction(INDEPENDENT_COST)+integral
    cost_upper = F.from_float(float(total_cost.hi))

    mtree_path = directory.parents[1]/'mtree_refinement'/'certificate.json'
    mtree = json.loads(mtree_path.read_text())['M_tree']
    m_lower = F(int(mtree['lower_numerator']),int(mtree['denominator']))
    deficit = 2*m_lower-cost_upper
    require(deficit > F(1,10000), 'comparison point did not prove the requested deficit')

    # Normalized masses have signed differences summing to zero. Each Gram
    # entry's generator values lie in [1/16,16], so the entry change is at
    # most ((16-1/16)/2)*||m-m_tree||_1 < 8*||m-m_tree||_1.
    mass_radius = F(1,10000)
    epsilon = 8*mass_radius
    require(epsilon < 1, 'first normalization not controlled')
    vertex_error = 2*epsilon/(1-epsilon)
    edge_error = 4*vertex_error
    edge_diagonal_lower = F(97,25)  # 2+2rho > 3.88, sqrt(2)>1.41.
    require(F(141,100)**2 < 2, 'square-root lower bracket failed')
    require(edge_error < edge_diagonal_lower, 'second normalization not controlled')
    correlation_error = 2*edge_error/(edge_diagonal_lower-edge_error)
    require(correlation_error < radius, 'full-simplex mass ball does not enter comparison region')

    # Every mixed Boolean difference is nonpositive: the Gaussian
    # monotonicity theorem applies along the PSD segment to the comparison.
    fourier, independent, max_difference = coefficients()
    require(independent==INDEPENDENT_COST, 'independent Boolean value changed')
    source_paths = [Path(__file__), directory/'p4_interval.py', arithmetic_path,
                    mtree_path]
    record = {
        'status':'RIGOROUS WIDE LOCAL NEIGHBORHOOD; NOT A WHOLE-CONE BOUND',
        'comparison_point':'exact natural P4 boundary covariance minus r in every off-diagonal entry',
        'radius':str(radius),
        'comparison_covariance_lo':[[float(x.lo) for x in row]for row in comparison],
        'comparison_covariance_hi':[[float(x.hi) for x in row]for row in comparison],
        'cholesky_pivots':pivots,
        'full_two_sided_box_shifted_cholesky_pivots':full_box_pivots,
        'total_cost_enclosure':exact_endpoints(total_cost),
        'total_cost_binary64_diagnostics':[float(total_cost.lo),float(total_cost.hi)],
        'total_cost_upper_fraction':str(cost_upper),
        'two_M_tree_minus_cost_exact_lower_bound':str(deficit),
        'two_M_tree_minus_cost_simple_strict_lower_bound':'1/10000',
        'one_sided_neighborhood':'all PSD unit-diagonal C with C_ij >= (C_tree)_ij-radius for every i<j',
        'full_simplex_L1_radius':str(mass_radius),
        'mass_ball_Gram_entry_error_bound':str(epsilon),
        'mass_ball_vertex_correlation_error_bound':str(vertex_error),
        'mass_ball_edge_Gram_entry_error_bound':str(edge_error),
        'mass_ball_key_correlation_error_bound':str(correlation_error),
        'integration':stats,
        'taylor_degree':args.degree,
        'rational_natural_masses':masses,
        'max_boolean_mixed_difference_absolute_value':max_difference,
        'upstream_interval_arithmetic_blob':arithmetic_blob,
        'input_sha256':{str(p.relative_to(repository)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in source_paths},
        'scope':[
            'Certified at the comparison point by the existing interval Taylor integrator.',
            'Extends to the stated one-sided region by Gaussian comparison for submodular costs.',
            'Contains the full mass-simplex L1 ball with radius 1/10000 around the natural point.',
            'Does not assert that every cone covariance lies in this region.'
        ]
    }
    text = json.dumps(record,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()
