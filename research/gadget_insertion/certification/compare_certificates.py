#!/usr/bin/env python3
"""Derive strict rational comparisons from already generated enclosures."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path


def require(condition,message):
    if not condition:raise ArithmeticError(message)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    directory=Path(__file__).resolve().parent
    candidate=json.loads((directory/'p4_point.json').read_text())
    tree=json.loads((directory/'p4_tree.json').read_text())
    mtree=json.loads((directory.parents[1]/'mtree_refinement'/'certificate.json').read_text())
    cl,cu=map(F,candidate['per_vertex_exact_endpoints']['exact_binary64_endpoint_fractions'])
    tl,tu=map(F,tree['per_vertex_exact_endpoints']['exact_binary64_endpoint_fractions'])
    ml=F(int(mtree['M_tree']['lower_numerator']),int(mtree['M_tree']['denominator']))
    mu=F(int(mtree['M_tree']['upper_numerator']),int(mtree['M_tree']['denominator']))
    gap=cl-tu
    require(gap>F(9,10**6),'candidate excess of 9e-6 per vertex was not certified')
    saving=2*ml-4*tu
    require(saving>F(433,50000),'natural total saving of 0.00866 was not certified')
    require(cu<F(6807,50000),'candidate point bound 0.13614 failed')
    for item in (candidate,tree):
        require(F(item['local_box_upper_endpoint_fraction'])<F(6807,50000),
                'local box bound 0.13614 failed')
        require(item['integration']['accepted_area_exact']=='1','incomplete integration cover')
    result={
        'comparisons_passed':True,
        'P4_natural_tree_is_not_a_maximizer_of_the_relaxed_cone':True,
        'candidate_minus_natural_per_vertex_exact_lower_bound':str(gap),
        'candidate_minus_natural_per_vertex_simple_strict_lower_bound':'9/1000000',
        'two_M_tree_minus_natural_P4_total_exact_lower_bound':str(saving),
        'two_M_tree_minus_natural_P4_total_simple_strict_lower_bound':'433/50000',
        'candidate_local_box_per_vertex_strict_upper_bound':'6807/50000',
        'candidate_local_correlation_radius':candidate['local_correlation_box_radius'],
        'candidate_local_simplex_L1_radius':candidate['local_simplex_L1_radius'],
        'scope':'Strict comparisons and local bounds only; no bound over the entire P4 cone.'}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':main()
