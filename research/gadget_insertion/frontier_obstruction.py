#!/usr/bin/env python3
"""Exact comparison between the reported frontier and the current P4 method.

No Gaussian integral is evaluated here. The lower cost endpoint is an input
from the separately replayable rigorous P4 point certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from path_cone import enumerate_types


def require(condition,message):
    if not condition:
        raise ArithmeticError(message)


def decimal(value,places=30,upper=False):
    scaled=value.numerator*10**places
    n=-((-scaled)//value.denominator) if upper else scaled//value.denominator
    return f'{n//10**places}.{n%10**places:0{places}d}'


def verify():
    root=Path(__file__).resolve().parent
    point_path=root/'certification'/'p4_point.json'
    tree_path=root/'certification'/'p4_tree.json'
    point=json.loads(point_path.read_text())
    tree=json.loads(tree_path.read_text())
    require(point['point']=='rationalized_nonuniform_candidate','wrong point certificate')
    require(tree['point']=='natural_tree','wrong natural certificate')
    type_set=set(enumerate_types(4))
    masses=point['rational_masses']
    require(all(tuple(row['type']) in type_set and F(row['mass'])>0 for row in masses),
            'candidate support is outside the declared P4 cone')
    require(sum((F(row['mass']) for row in masses),F(0))==1,'candidate masses do not sum to one')
    lower,upper=map(F,point['per_vertex_exact_endpoints']['exact_binary64_endpoint_fractions'])
    natural_lower=F(tree['per_vertex_exact_endpoints']['exact_binary64_endpoint_fractions'][0])
    require(0<lower<=upper,'invalid candidate interval')

    reported_frontier=F(6341,1341)
    matching_kappa=1/(2*(reported_frontier-1))
    require(matching_kappa==F(1341,10000),'frontier conversion failed')
    require(lower>matching_kappa,'candidate does not obstruct the target')
    require(natural_lower>matching_kappa,'natural point does not obstruct the target')
    method_ceiling=1+1/(2*lower)
    require(method_ceiling<reported_frontier,'method ceiling does not lie below frontier')
    former_target=F(6807,50000)
    former_consequence=1+1/(2*former_target)
    require(former_consequence==F(31807,6807),'former target conversion failed')
    require(former_consequence<reported_frontier,'former target not superseded')

    return {
        'status':'EXACT METHOD-OBSTRUCTION COMPARISONS PASSED',
        'reported_frontier':{
            'source_url':'https://www.yukon.org/circuitbreaker',
            'retrieval_date':'2026-11-17',
            'reported_circuit_coefficient':str(reported_frontier),
            'decimal_lower':decimal(reported_frontier),
            'decimal_upper':decimal(reported_frontier,upper=True),
            'source_status':'Reported provisional result, pending full mathematical review; leading proof not independently verified here.',
        },
        'cutwidth_coefficient_to_match_reported_frontier':str(matching_kappa),
        'to_exceed_reported_frontier':'Need kappa < 1341/10000 under c=1+1/(2*kappa).',
        'certified_candidate_per_vertex_lower':str(lower),
        'certified_candidate_per_vertex_decimal_lower':decimal(lower),
        'certified_natural_per_vertex_lower':str(natural_lower),
        'candidate_excess_over_frontier_kappa':str(lower-matching_kappa),
        'uniform_P4_cone_method_circuit_ceiling':str(method_ceiling),
        'uniform_P4_cone_method_circuit_decimal_ceiling':decimal(method_ceiling,upper=True),
        'reported_frontier_minus_method_ceiling':str(reported_frontier-method_ceiling),
        'reported_frontier_minus_method_ceiling_decimal_lower':decimal(reported_frontier-method_ceiling),
        'former_target_kappa':str(former_target),
        'former_target_circuit_coefficient':str(former_consequence),
        'former_target_would_not_reach_reported_frontier':True,
        'natural_P4_cost_also_exceeds_frontier_target':True,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (point_path,tree_path)},
        'scope':[
            'An upper ceiling on coefficients obtained by a uniform bound over the present relaxed P4 cone followed by the stated reduction.',
            'Not an upper bound on circuit complexity or on all possible graph-ordering methods.',
            'The rational candidate need not be realizable as a finite graph neighborhood.',
            'Additional realizability or averaging arguments can change the proof obligation; changing the gadget or insertion objective can also evade this obstruction.',
            'The public frontier is reported, not independently mathematically certified by this comparison.'
        ]
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
