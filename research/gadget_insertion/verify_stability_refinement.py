#!/usr/bin/env python3
"""Exact finite/scalar checks for STABILITY_REFINEMENT.md.

These checks support the written analytic proof and its stated upstream and
wide-neighborhood certificate inputs. They do not replace those inputs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from path_cone import enumerate_types, structure
from path_gadgets import Qsqrt2 as Q
from verify_rigidity import exact_rref, power_two


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def verify():
    checks = []
    def check(name, condition):
        require(condition, name)
        checks.append(name)

    sqrt2_lo, sqrt2_hi = F(141421356237,10**11), F(141421356238,10**11)
    check('sqrt2_enclosure', sqrt2_lo**2 < 2 < sqrt2_hi**2)
    def lower(z):
        return z.a+z.b*(sqrt2_lo if z.b>=0 else sqrt2_hi)
    def power(exponent):
        coefficient=F(2)**(exponent//2)
        return Q(0,coefficient) if exponent%2 else Q(coefficient,0)

    # On the old LOC chart the squared cross-entry gradient norm is at most
    # this rational number. Comparing it to 16*lambda converts the old
    # quadratic cost deficit to a 4*sqrt(eta) covariance radius.
    chart_gradient_squared = 2*(F(5,4)*(1+F(1,7500)))**2+2*(1+F(1,10000))**2
    check('chart_gradient_norm', chart_gradient_squared < 16*F(3349,10000))
    check('LOC_S_bounds', 2+sqrt2_lo-F(1,1000) > F(9,4)
          and 2+sqrt2_hi+F(1,1000) < 4)

    eta=F(17,5*10**19)
    v_upper=F(729,1250000000000)
    r_upper=5*v_upper
    u_upper=F(27,500000)
    check('sqrt_eta_upper', eta < v_upper**2)
    check('exact_sqrt_r_upper', r_upper == u_upper**2)
    check('monotone_path_remains_local', 4*v_upper < F(1,2000))
    check('monotone_path_loss', eta < F(1,200000))
    check('pair_covariance_radius', 100*v_upper <= 1)

    u0=F(1,10**4)
    check('small_u_regime', u_upper < u0)
    check('row_norm_upper', (1+5*u0)**4 < F(1003,1000))
    check('squared_row_norm_upper', F(1003,1000)**2 < F(1007,1000))
    check('raw_residual_coefficient', F(1003,1000)*(21+8*u0) < 22)
    check('c0_lower_bound', (1+2*sqrt2_lo/3)/2 > F(98,100)**2)
    check('aggregate_edge_excess_constant', F(32,3)/F(98,100) < 11)
    check('aggregate_c_excess_constant', F(11)/(4*F(98,100)) < 3)
    check('aggregate_normalized_residual_energy', (8+12*u0)**2 < 65)
    check('zero_edge_EC_slack_constant', F(3)/sqrt2_hi-2 > F(3,25))
    check('zero_edge_mass_constant', 8*F(15,7)*F(1007,1000)*11/F(3,25) < 1600)
    check('nonzero_increment_ratio_separation', 1/sqrt2_hi > F(7,10))
    check('normalization_ratio_error', 15*sqrt2_hi < 22)
    check('bad_pm1_normalized_energy_floor',
          (F(7,10)-22*u0)**2/(8*F(1007,1000)) > F(3,50))
    check('bad_pm1_mass_constant', F(65)/F(3,50) < 1100)
    for depth in range(1,5):
        check(f'diagonal_error_depth_{depth}',
              (1+5*u0)**(2*depth)-1 <= F(1002,100)*depth*u0)

    types=enumerate_types(4)
    owners, edges, _, neighbors=structure(4)
    root2=Q(0,1)
    energy_floor=Q(F(3,16),F(-1,8))
    tree_indices=[]
    floor_indices=[]
    for index,t in enumerate(types):
        g=[power(-value) for value in t]
        energy=Q()
        for v in range(4):
            residual=sum((g[w] for w in neighbors[v]),Q())+(-2)*root2*g[v]
            energy+=residual*residual
        if energy.a==energy.b==0:
            tree_indices.append(index)
        else:
            difference=energy+(-1)*energy_floor
            if difference.a==difference.b==0:
                floor_indices.append(index)
            else:
                require(lower(difference)>0,'non-tree energy floor failed')
    check('tree_zero_set', tree_indices==[455,457,904,1794,4438,4600])
    check('core_generator_square_minimum', all(t[v]<=3 for t in types for v in range(4)))
    check('positive_energy_minimum_attained', bool(floor_indices))
    check('off_face_mass_constant', lower(energy_floor) > F(4*22**2,200000))

    rows=(0,1,2,3,4,8)
    depths=(0,1,2,3,1,4)
    matrix=[[power_two(-types[index][row]) for index in tree_indices] for row in rows]
    reduced,pivots=exact_rref([row+[F(i==j) for j in range(6)]
                              for i,row in enumerate(matrix)])
    require(pivots==list(range(6)),'tree minor not invertible')
    inverse=[row[6:] for row in reduced]
    column_sums=[sum(abs(inverse[i][j]) for i in range(6)) for j in range(6)]
    weighted_sum=sum(c*d for c,d in zip(column_sums,depths))
    check('depth_weighted_inverse_norm', weighted_sum==F(79,12))
    check('diagonal_contribution_constant', F(1002,100)*weighted_sum < 66)
    off_face_coefficient=F(0)
    off_face_maximizers=[]
    for index,t in enumerate(types):
        if index in tree_indices:
            continue
        diagonal=[power_two(-t[row]) for row in rows]
        mapped=[sum(a*b for a,b in zip(row,diagonal)) for row in inverse]
        coefficient=1+sum(map(abs,mapped))
        if coefficient>off_face_coefficient:
            off_face_coefficient=coefficient
            off_face_maximizers=[index]
        elif coefficient==off_face_coefficient:
            off_face_maximizers.append(index)
    check('off_face_recovery_coefficient', off_face_coefficient==F(9,2))

    mass_upper=66*u_upper+12150*r_upper
    check('inside_wide_mass_neighborhood', mass_upper < F(9,2500))
    root=Path(__file__).resolve().parent
    wide_path=root/'certification'/'wide_neighborhood.json'
    mass_path=root/'certification'/'mass_neighborhood.json'
    tree_path=root.parent/'mtree_refinement'/'certificate.json'
    wide=json.loads(wide_path.read_text())
    mass_certificate=json.loads(mass_path.read_text())
    tree=json.loads(tree_path.read_text())['M_tree']
    tree_lower=F(int(tree['lower_numerator']),int(tree['denominator']))
    check('wide_mass_radius_input', F(mass_certificate['full_simplex_L1_radius']) >= F(9,2500))
    check('wide_mass_correlation_implication',
          F(mass_certificate['maximum_downward_correlation_loss_fraction'])
          < F(wide['radius']) == F(1,300))
    check('wide_mass_cost_input', F(mass_certificate['total_P4_cost_upper_fraction'])
          == F(wide['total_cost_upper_fraction']))
    check('wide_deficit_input',
          2*tree_lower-F(wide['total_cost_upper_fraction']) > F(1,10000) > eta)
    check('upstream_nonlocal_gap',
          tree_lower-(F('0.27658360687864314')+F(1,10**13)) > F(3,10**12) > eta)

    return {
        'status':'EXACT FINITE AND SCALAR CHECKS PASSED; SEE ANALYTIC PROOF',
        'checks':checks,
        'global_P4_deficit':str(eta),
        'pair_covariance_radius_rule':'5*sqrt(eta)',
        'pair_covariance_radius_upper':str(r_upper),
        'sqrt_pair_radius_upper':str(u_upper),
        'full_simplex_mass_distance_upper':str(mass_upper),
        'required_mass_neighborhood_radius':'9/2500',
        'tree_support_indices':tree_indices,
        'off_face_minimum_energy':energy_floor.json(),
        'off_face_minimum_energy_attained_at':floor_indices,
        'raw_energy_off_face_mass_bound':'200000*r',
        'refined_off_face_mass_bound':'2700*r',
        'tree_diagonal_minor_rows':list(rows),
        'tree_diagonal_minor_inverse':[[str(x) for x in row] for row in inverse],
        'inverse_absolute_column_sums':[str(x) for x in column_sums],
        'depth_weighted_column_sum':str(weighted_sum),
        'off_face_recovery_coefficient':str(off_face_coefficient),
        'off_face_recovery_maximizers':off_face_maximizers,
        'mass_distance_rule':'66*sqrt(r)+12150*r',
        'input_sha256':{str(path.relative_to(root.parent)):hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in (wide_path,mass_path,tree_path)},
        'upstream_finite_certificates_replayed':False,
        'new_leaderboard_record':False,
        'scope':'Supports STABILITY_REFINEMENT.md; requires its upstream analytic inputs and the separately verified wide natural neighborhood.'
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=verify()
    text=json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()
