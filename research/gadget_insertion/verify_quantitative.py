#!/usr/bin/env python3
"""Exact arithmetic checks for the conservative P4 global-gap proof.

This checks the finite linear algebra and scalar implications used in
QUANTITATIVE_GLOBAL.md. It is not a replacement for that analytic proof or for
the upstream pair-theorem certificates on which the proof depends.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import sys

from path_cone import enumerate_types, structure
from verify_rigidity import exact_rref, power_two


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def decimal_bound(value, places=60, upper=False):
    scaled = value.numerator*10**places
    n = -((-scaled)//value.denominator) if upper else scaled//value.denominator
    return f'{n//10**places}.{n%10**places:0{places}d}'


def verify():
    root = Path(__file__).resolve().parent
    checks = []

    def check(name, condition):
        require(condition, name)
        checks.append(name)

    # Every comparison below is rational. The square checks establish the
    # stated bracket for sqrt(2), independently of floating-point arithmetic.
    sqrt2_lo = Q(141421356237, 10**11)
    sqrt2_hi = Q(141421356238, 10**11)
    check('sqrt2_bracket', sqrt2_lo**2 < 2 < sqrt2_hi**2)
    p_hi = (1 + 2*sqrt2_hi)/4
    check('old_tree_box_minimum_eigenvalue', 1-p_hi-Q(3,1000) > Q(3,100))
    check('pi_upper_square', Q(22,7)**2 < 10)
    check('cross_derivative_magnitude_margin', Q(1,5)**2 < Q(6,100))

    eta = Q(1, 10**54)
    sqrt_eta = Q(1, 10**27)
    r = 32*sqrt_eta
    z_upper = Q(2, 10**13)
    check('exact_sqrt_eta', sqrt_eta**2 == eta)
    check('sqrt_r_upper', r < z_upper**2)
    check('small_r_regime', r <= Q(1, 10**20))

    mtree_path = root.parent/'mtree_refinement'/'certificate.json'
    mtree = json.loads(mtree_path.read_text())['M_tree']
    mtree_lo = Q(int(mtree['lower_numerator']), int(mtree['denominator']))
    mtree_hi = Q(int(mtree['upper_numerator']), int(mtree['denominator']))
    old_lower = Q('0.27658360687864314')
    check('outside_LOC_gap', mtree_lo-(old_lower+Q(1,10**13)) > Q(3, 10**12) > eta)
    check('LOC_norm_implication', 4*Q('0.3349') > 1)
    check('tight_cut_C0_within_half_box', 16*sqrt_eta < Q(1,2000))
    check('monotone_box_exit_loss_exceeds_eta',
          Q(1,100)*(Q(1,1000)-16*sqrt_eta) > eta)
    check('pair_covariance_bound', 16*sqrt_eta+100*eta <= r)

    # Edge norm ratios satisfy |u-1|<=6 sqrt(r). All extended-path
    # vertices are at most four required edges from root zero.
    z_max = Q(1,10**10)
    check('edge_ratio_constant', 32 < 6**2)
    check('four_edge_norm_propagation', (1+6*z_max)**4-1 <= 30*z_max)
    check('squared_diagonal_bound', 60+900*z_max <= 61)
    check('raw_star_residual_bound', 8+90+60*sqrt2_hi < 200)

    types = enumerate_types(4)
    owners, edges, _, neighbors = structure(4)
    check('cone_size', len(types) == 4655)
    check('generator_range', all(max(map(abs,t)) <= 4 for t in types))
    equality_indices = []
    residual_patterns = set()
    for index,t in enumerate(types):
        good = True
        for v in range(4):
            increments = tuple(t[u]-t[v] for u in neighbors[v])
            require(all(d in (-1,0,1) for d in increments), 'invalid edge increment')
            a = Q(increments.count(0))
            b = Q(increments.count(-1))+Q(increments.count(1),2)-2
            if a == b == 0:
                require(sorted(increments) == [-1,1,1], 'unexpected zero star residual')
            else:
                good = False
                lower = a+b*(sqrt2_lo if b>=0 else sqrt2_hi)
                upper = a+b*(sqrt2_hi if b>=0 else sqrt2_lo)
                require(lower >= Q(1,10) or upper <= -Q(1,10),
                        'star residual below the proved separation')
                residual_patterns.add(increments)
        if good:
            equality_indices.append(index)
    check('only_six_zero_residual_types', len(equality_indices)==6)

    selected_rows = (0,1,2,3,4,8)
    tree_types = [types[i] for i in equality_indices]
    matrix = [[power_two(-t[i]) for t in tree_types] for i in selected_rows]
    augmented = [row+[Q(i==j) for j in range(6)] for i,row in enumerate(matrix)]
    reduced, pivots = exact_rref(augmented)
    require(pivots == list(range(6)), 'tree diagonal minor is singular')
    inverse = [row[6:] for row in reduced]
    for i,j in itertools.product(range(6), repeat=2):
        require(sum(matrix[i][k]*inverse[k][j] for k in range(6)) == int(i==j),
                'matrix inverse check failed')
    inverse_norm = max(sum(abs(x) for x in row) for row in inverse)
    check('tree_minor_inverse_norm', inverse_norm == 3)

    check('off_face_mass_union_bound', 4*200**2*40**2 <= 3*10**8)
    check('mass_distance_bound', 1098+289*(3*10**8)*z_max <= 2000)
    mass_distance_upper = 2000*z_upper
    check('inside_natural_mass_neighborhood', mass_distance_upper < Q(1,2**30))

    natural_path = root/'certification'/'p4_tree.json'
    natural = json.loads(natural_path.read_text())
    check('natural_certificate_identity', natural['point']=='natural_tree')
    check('natural_mass_radius', Q(natural['local_simplex_L1_radius'])==Q(1,2**30))
    local_upper = Q(natural['local_box_upper_endpoint_fraction'])
    target = Q(6807,50000)
    check('natural_local_target', local_upper < target)
    check('local_branch_strict_global_gap', 4*target < 2*mtree_lo-eta)
    check('positive_new_cutwidth_bound', mtree_lo/2-eta/7 > 0)

    # Re-run the independent integer certificate at sufficient precision to
    # exhibit the extremely small improvement as explicit rational bounds.
    # Its default 30-place enclosure is insufficient for this comparison.
    sys.path.insert(0,str(root.parent/'mtree_refinement'))
    from certify import certify
    refined = certify(digits=120,terms=640,output_digits=60)
    high = refined['M_tree']
    m_lo = Q(int(high['lower_numerator']),int(high['denominator']))
    m_hi = Q(int(high['upper_numerator']),int(high['denominator']))
    kappa_upper = m_hi/2-eta/7
    circuit_lower = 1+1/(m_hi-2*eta/7)
    old_circuit_upper = 1+1/m_lo
    check('new_rational_cutwidth_bound_below_old_exact_value', kappa_upper < m_lo/2)
    check('new_circuit_coefficient_gain_exceeds_3eta', circuit_lower-old_circuit_upper > 3*eta)

    return {
        'status':'EXACT FINITE AND SCALAR CHECKS PASSED; ANALYTIC PROOF REQUIRED',
        'checks':checks,
        'delta_lower_bound':str(eta),
        'pair_covariance_radius':str(r),
        'mass_distance_upper_bound':str(mass_distance_upper),
        'natural_mass_neighborhood_radius':natural['local_simplex_L1_radius'],
        'tree_support_indices':equality_indices,
        'tree_diagonal_minor_rows':list(selected_rows),
        'tree_diagonal_minor_inverse':[[str(x) for x in row] for row in inverse],
        'tree_diagonal_minor_inverse_infinity_norm':str(inverse_norm),
        'nonzero_residual_patterns_checked':len(residual_patterns),
        'input_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (mtree_path,natural_path)},
        'refined_tree_integer_certificate':refined,
        'consequence_of_the_analytic_global_proof':{
            'cubic_cutwidth_coefficient_upper_fraction':str(kappa_upper),
            'cubic_cutwidth_coefficient_decimal_upper':decimal_bound(kappa_upper,upper=True),
            'circuit_coefficient_lower_fraction':str(circuit_lower),
            'circuit_coefficient_decimal_lower':decimal_bound(circuit_lower),
            'strict_gain_over_original_exact_circuit_coefficient':str(3*eta),
            'new_leaderboard_record':False,
            'full_upstream_certificates_replayed':False,
        },
        'scope':[
            'Exact finite checks supporting QUANTITATIVE_GLOBAL.md, not a standalone proof.',
            'Uses upstream LOC/R1 and pair-theorem certificates as stated dependencies.',
            'Uses the separately verified natural P4 local certificate.',
            'The resulting global gap is tiny and establishes no leaderboard record.'
        ]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
