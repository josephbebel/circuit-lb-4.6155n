#!/usr/bin/env python3
"""Exact numerical consequence of the improved stability and P4-only proof.

The analytic graph and Gaussian proofs remain explicit dependencies. This
script replays the small exact stability checks and the integer tree-value
enclosure, then translates their proved inequalities without float rounding.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

from verify_stability_refinement import verify


def require(condition,message):
    if not condition:
        raise ArithmeticError(message)


def decimal_bound(value,places=45,upper=False):
    scaled=value.numerator*10**places
    n=-((-scaled)//value.denominator) if upper else scaled//value.denominator
    return f'{n//10**places}.{n%10**places:0{places}d}'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    stability=verify()
    eta=F(stability['global_P4_deficit'])
    sys.path.insert(0,str(root.parent/'mtree_refinement'))
    from certify import certify
    tree=certify(digits=90,terms=480,output_digits=45)
    bounds=tree['M_tree']
    lower=F(int(bounds['lower_numerator']),int(bounds['denominator']))
    upper=F(int(bounds['upper_numerator']),int(bounds['denominator']))
    # GLOBAL_REDUCTION.md proves kappa <= M4/4, and the new stability
    # theorem proves M4 <= 2*M_tree-eta.
    p4_upper=2*upper-eta
    kappa_upper=p4_upper/4
    circuit_lower=1+1/(2*kappa_upper)
    old_circuit_upper=1+1/lower
    require(kappa_upper>0,'nonpositive cutwidth coefficient')
    require(kappa_upper<lower/2,'new bound does not improve the exact old method')
    require(circuit_lower-old_circuit_upper>F(22,10**19),
            'claimed strict circuit-coefficient improvement not verified')
    paths=[root/'STABILITY_REFINEMENT.md',root/'GLOBAL_REDUCTION.md',
           root/'verify_stability_refinement.py',
           root/'certification'/'mass_neighborhood.json',
           root/'certification'/'wide_neighborhood.json',
           root.parent/'mtree_refinement'/'certify.py']
    record={
        'status':'EXACT CONSEQUENCE OF THE STATED ANALYTIC PROOFS',
        'global_P4_deficit':str(eta),
        'stability_scalar_and_finite_checks':len(stability['checks']),
        'global_cutwidth_reduction':'kappa <= M4/4',
        'M4_total_cost_upper_fraction':str(p4_upper),
        'M4_total_cost_decimal_upper':decimal_bound(p4_upper,upper=True),
        'cubic_cutwidth_coefficient_upper_fraction':str(kappa_upper),
        'cubic_cutwidth_coefficient_decimal_upper':decimal_bound(kappa_upper,upper=True),
        'circuit_coefficient_lower_fraction':str(circuit_lower),
        'circuit_coefficient_decimal_lower':decimal_bound(circuit_lower),
        'cutwidth_decrease_from_old_exact_method':str(eta/4),
        'circuit_gain_over_old_exact_method_strict_lower_bound':str(F(22,10**19)),
        'tree_integer_certificate':tree,
        'input_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in paths},
        'scope':[
            'Uniform asymptotic kernel cutwidth coefficient and upstream all-lengths circuit consequence.',
            'Requires STABILITY_REFINEMENT.md, GLOBAL_REDUCTION.md and their named upstream results.',
            'The small new certificates have replay scripts; the large upstream suites were not replayed.',
            'This does not reach the updated reported frontier or establish a new leaderboard record.'
        ],
        'new_leaderboard_record':False,
        'full_upstream_certificates_replayed':False,
    }
    text=json.dumps(record,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()
