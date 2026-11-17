#!/usr/bin/env python3
"""Interval first-derivative proof of a natural-P4 mass neighborhood.

This bounds the entire normalized-simplex L1 ball. It invokes the separately
certified one-sided covariance comparison; it performs no quadrature.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from p4_interval import IA, interval_fraction, require


def minimum_interval(a, axis=0):
    return IA(np.min(a.lo,axis=axis),np.min(a.hi,axis=axis))


def maximum_interval(a, axis=0):
    return IA(np.max(a.lo,axis=axis),np.max(a.hi,axis=axis))


def clipped_correlation(value):
    # All actual Gram entries are positive and the matrices are PSD.
    return IA(np.maximum(value.lo,0),np.minimum(value.hi,1))


def certify(radius):
    here=Path(__file__).resolve().parent
    radius=F(radius)
    require(0<radius<=F(1,100),'radius outside checked small-ball range')
    arithmetic_path=here.parents[2]/'mtree'/'enclosure'/'ivarith.py'
    arithmetic=arithmetic_path.read_bytes()
    arithmetic_blob=hashlib.sha1(b'blob '+str(len(arithmetic)).encode()+b'\0'+arithmetic).hexdigest()
    require(arithmetic_blob=='b99d01f757f2b3d75151ba9403a547c6b8c60477',
            'interval arithmetic differs from the audited upstream version')
    sys.path.insert(0,str(here.parent))
    from path_cone import enumerate_types, tree_distribution, boundary_owners
    types=enumerate_types(4)
    require(len(types)==4655,'incorrect cone size')
    t=np.asarray(types,dtype=np.int64)
    require(bool(np.all(t[:,0]==0) and np.max(np.abs(t))<=4),'generator range failed')
    masses=tree_distribution(4,types)
    support=np.flatnonzero([bool(m) for m in masses])
    exponent=-t[:,:,None]-t[:,None,:]
    # ldexp is exact for these small integral exponents; no transcendental
    # exp2 implementation is used as an exact constant source.
    coefficient=IA(np.ldexp(np.ones_like(exponent,dtype=float),exponent//2))*IA(
        np.where(exponent%2,float(IA(2).sqrt().lo),1.0),
        np.where(exponent%2,float(IA(2).sqrt().hi),1.0))
    center=IA(np.zeros((10,10)))
    for i in support:
        center+=coefficient[i]*interval_fraction(masses[i])
    # Exact natural Gram diagonals equal one. Preserve this identity in
    # the center enclosure rather than accumulating artificial width.
    for j in range(10):
        require(sum(masses[i]*F(2)**(-types[i][j]) for i in support)==1,
                'natural diagonal identity failed')
        center.lo[j,j]=center.hi[j,j]=1

    # On the path from m_tree to m, negative mass changes occur only at
    # tree support indices. Each signed total is ||m-m_tree||_1/2.
    all_min=minimum_interval(coefficient)
    all_max=maximum_interval(coefficient)
    support_min=minimum_interval(coefficient[support])
    support_max=maximum_interval(coefficient[support])
    spread=IA(all_min.lo,all_max.hi)-IA(support_min.lo,support_max.hi)
    gram=center+interval_fraction(radius/2)*spread
    require(bool(np.all(np.diag(gram.lo)>0)),'Gram diagonal box is not positive')
    diagonal=[gram[i,i] for i in range(10)]
    vertex={}
    for i in range(10):
        for j in range(i,10):
            vertex[i,j]=IA(1) if i==j else clipped_correlation(
                gram[i,j]/(diagonal[i]*diagonal[j]).sqrt())
    direction=coefficient-center
    derivative={}
    for i in range(10):
        for j in range(i,10):
            derivative[i,j]=IA(np.zeros(len(types))) if i==j else (
                direction[:,i,j]/(diagonal[i]*diagonal[j]).sqrt()
                -vertex[i,j]*(direction[:,i,i]/diagonal[i]
                               +direction[:,j,j]/diagonal[j])/2)
    def at(table,i,j):
        return table[min(i,j),max(i,j)]
    edges=[(owner,4+j) for j,owner in enumerate(boundary_owners(4))]
    edge_diagonal=[2+2*at(vertex,*e) for e in edges]
    edge_derivative=[2*at(derivative,*e) for e in edges]
    outcomes=[]
    maximum_loss=F(0)
    for i in range(6):
        for j in range(i+1,6):
            w=sum((at(vertex,a,b) for a in edges[i] for b in edges[j]),IA(0))
            dw=sum((at(derivative,a,b) for a in edges[i] for b in edges[j]),IA(0))
            denominator=(edge_diagonal[i]*edge_diagonal[j]).sqrt()
            corr=clipped_correlation(w/denominator)
            gradients=dw/denominator-corr*(edge_derivative[i]/edge_diagonal[i]
                                         +edge_derivative[j]/edge_diagonal[j])/2
            # Along the segment, dC= sum_t gradient_t * (m_t-mT_t).
            # The negative part is confined to the original tree support.
            slope=(maximum_interval(gradients[support])-minimum_interval(gradients))/2
            slope_upper=F.from_float(float(slope.hi))
            loss=radius*max(F(0),slope_upper)
            maximum_loss=max(maximum_loss,loss)
            outcomes.append({'pair':[i,j],
                             'downward_slope_upper_fraction':str(slope_upper),
                             'correlation_loss_upper_fraction':str(loss),
                             'slope_binary64_diagnostic':float(slope.hi)})
    wide_path=here/'wide_neighborhood.json'
    wide=json.loads(wide_path.read_text())
    comparison_radius=F(wide['radius'])
    require(comparison_radius==F(1,300),'unexpected comparison certificate radius')
    mtree_path=here.parents[1]/'mtree_refinement'/'certificate.json'
    mtree=json.loads(mtree_path.read_text())['M_tree']
    m_lower=F(int(mtree['lower_numerator']),int(mtree['denominator']))
    require(2*m_lower-F(wide['total_cost_upper_fraction'])>F(1,10000),
            'wide certificate does not have the required exact deficit')
    passed=maximum_loss<comparison_radius
    paths=[Path(__file__),here/'p4_interval.py',wide_path,here.parent/'path_cone.py',mtree_path]
    return {
        'status':'RIGOROUS MASS-NEIGHBORHOOD DERIVATIVE CHECK' if passed else 'REGION NOT CERTIFIED',
        'certified':passed,
        'full_simplex_L1_radius':str(radius),
        'comparison_correlation_radius':str(comparison_radius),
        'maximum_downward_correlation_loss_fraction':str(maximum_loss),
        'maximum_downward_correlation_loss_diagnostic':float(maximum_loss),
        'pairs':outcomes,
        'total_P4_cost_upper_fraction':wide['total_cost_upper_fraction'],
        'two_M_tree_minus_cost_simple_strict_lower_bound':wide['two_M_tree_minus_cost_simple_strict_lower_bound'],
        'source_types':len(types),
        'tree_support_indices':support.tolist(),
        'upstream_interval_arithmetic_blob':arithmetic_blob,
        'input_sha256':{str(p.relative_to(here.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in paths},
        'scope':'Full normalized mass simplex L1 ball about exact natural masses; not the entire cone.',
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',default='9/2500')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    radius=F(args.radius)
    require(0<radius<=F(1,100),'radius outside checked small-ball range')
    result=certify(radius)
    text=json.dumps(result,indent=2)+'\n'
    if args.output:
        require(result['certified'],'requested region was not certified')
        args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()
