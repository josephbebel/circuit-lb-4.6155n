#!/usr/bin/env python3
"""Exact Boolean cycle/path identities and finite alternating-window bounds.

No Gaussian integration or universal cutwidth claim is performed here.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path


def require(condition,message):
    if not condition:
        raise ArithmeticError(message)


def change_edges(colors,cycle=False):
    n=len(colors)
    return [int(colors[i]!=colors[(i+1)%n]) for i in range(n if cycle else n-1)]


def path_matching(changes):
    previous=0
    total=0
    for value in changes:
        previous=value*(1-previous)
        total+=previous
    return total


def cycle_matching(changes):
    n=len(changes)
    if all(changes):
        require(n%2==0,'an odd cycle cannot change sign on every edge')
        return n//2
    cut=changes.index(0)
    return path_matching([changes[(cut+1+i)%n] for i in range(n)])


def window_observable(colors):
    """H_K on K+1 ordered colors, using windows ending at the final edge."""
    changes=change_edges(colors)
    product=1
    total=0
    for length,value in enumerate(reversed(changes),1):
        product*=value
        total+=(1 if length%2 else -1)*product
    return total


def cycle_window_sum(colors,order):
    n=len(colors)
    return sum(window_observable([colors[(i-order+j)%n] for j in range(order+1)])
               for i in range(n))


def cycle_alternating_count(colors,order):
    changes=change_edges(colors,True)
    n=len(changes)
    return sum(all(changes[(i-j)%n] for j in range(order)) for i in range(n))


def brute_split_cost(colors,cycle=False):
    """One fixed external bit per vertex; path endpoints have one free edge."""
    n=len(colors)
    edges=n if cycle else n-1
    result=n
    for free in itertools.product((0,1),repeat=edges):
        value=0
        for i,color in enumerate(colors):
            incident=[color]
            if cycle or i:
                incident.append(free[(i-1)%edges])
            if cycle or i<n-1:
                incident.append(free[i%edges])
            value+=int(min(incident)!=max(incident))
        result=min(result,value)
    return result


def cubic_path_cost(boundary):
    """Actual P_k: two boundary bits at each endpoint, one at each interior."""
    k=len(boundary)-2
    require(k>=2,'path must contain at least two vertices')
    colors=[boundary[0]]+list(boundary[2:k])+[boundary[k]]
    forced=[False]*k
    forced[0]=boundary[0]!=boundary[1]
    forced[-1]=boundary[-2]!=boundary[-1]
    changes=[int(not forced[i] and not forced[i+1] and colors[i]!=colors[i+1])
             for i in range(k-1)]
    return sum(forced)+path_matching(changes)


def cubic_path_window_upper(boundary,order):
    k=len(boundary)-2
    require(order>=1 and order%2==1 and k>=order+1,'invalid path window order')
    colors=[boundary[0]]+list(boundary[2:k])+[boundary[k]]
    return (order+3)//2+sum(window_observable(colors[j:j+order+1])
                            for j in range(k-order))


def walsh_numerators(table):
    result=list(table)
    stride=1
    while stride<len(result):
        for start in range(0,len(result),2*stride):
            for i in range(stride):
                a,b=result[start+i],result[start+stride+i]
                result[start+i],result[start+stride+i]=a+b,a-b
        stride*=2
    return result


def observable_record(order):
    n=order+1
    table=[window_observable([(mask>>i)&1 for i in range(n)]) for mask in range(1<<n)]
    require(set(table)<={0,1},'window observable not Boolean')
    require(all(table[m]==table[m^((1<<n)-1)] for m in range(1<<n)),
            'window observable lacks complement symmetry')
    coefficients=walsh_numerators(table)
    require(all(not coefficients[m] for m in range(1<<n) if m.bit_count()%2),
            'odd Walsh coefficient survived')
    minimum,maximum=0,0
    for i,j in itertools.combinations(range(n),2):
        for mask in range(1<<n):
            if mask&((1<<i)|(1<<j)):
                continue
            delta=(table[mask|(1<<i)|(1<<j)]-table[mask|(1<<i)]
                   -table[mask|(1<<j)]+table[mask])
            minimum=min(minimum,delta);maximum=max(maximum,delta)
    return {
        'order':order,'boundary_variables':n,'patterns':len(table),
        'table_little_endian':table,
        'Walsh_denominator':len(table),
        'nonzero_Walsh_coefficients':[
            {'mask':mask,'numerator':value,'degree':mask.bit_count()}
            for mask,value in enumerate(coefficients) if value],
        'maximum_Walsh_degree':max(mask.bit_count() for mask,value in enumerate(coefficients) if value),
        'minimum_mixed_difference':minimum,'maximum_mixed_difference':maximum,
        'submodular':maximum<=0,
    }


def verify(max_cycle_size=8,max_path_size=8):
    from path_gadgets import cost_table
    path_checks=cycle_checks=window_checks=cubic_path_checks=endpoint_window_checks=0
    for n in range(2,max_cycle_size+1):
        for colors in itertools.product((0,1),repeat=n):
            changes=change_edges(colors,True)
            expected=cycle_matching(changes)
            require(brute_split_cost(colors,True)==expected,'cycle vertex-cover identity failed')
            cycle_checks+=1
            for order in range(1,10):
                value=cycle_window_sum(colors,order)
                if order%2:
                    require(expected<=value,'odd window upper bound failed')
                    require(value-expected<=cycle_alternating_count(colors,order),
                            'odd window error bound failed')
                else:
                    require(value<=expected,'even window lower bound failed')
                window_checks+=1
    for n in range(2,max_path_size+1):
        for colors in itertools.product((0,1),repeat=n):
            require(brute_split_cost(colors)==path_matching(change_edges(colors)),
                    'path vertex-cover identity failed')
            path_checks+=1
        table=cost_table(n)
        for mask in range(1<<(n+2)):
            boundary=[(mask>>i)&1 for i in range(n+2)]
            require(cubic_path_cost(boundary)==int(table[mask]),'cubic endpoint identity failed')
            cubic_path_checks+=1
            for order in range(1,min(n,10),2):
                require(cubic_path_window_upper(boundary,order)>=int(table[mask]),
                        'bounded path endpoint/window inequality failed')
                endpoint_window_checks+=1
    return {
        'status':'EXACT BOOLEAN IDENTITIES AND WINDOW BOUNDS PASSED',
        'one_boundary_bit_path_patterns_checked':path_checks,
        'one_boundary_bit_cycle_patterns_checked':cycle_checks,
        'cyclic_window_inequalities_checked':window_checks,
        'cubic_path_boundary_patterns_checked':cubic_path_checks,
        'bounded_path_window_inequalities_checked':endpoint_window_checks,
        'max_bruteforce_cycle_vertices':max_cycle_size,
        'max_bruteforce_path_vertices':max_path_size,
        'window_observables':[observable_record(k) for k in (1,3,5,7,9)],
        'scope':[
            'Exact Boolean identities and finite-window inequalities only.',
            'Odd window truncations are pointwise upper bounds; no independence is assumed.',
            'The individual odd window observable is generally not submodular.',
            'No Gaussian expectation, whole-cone bound, concentration theorem, or new circuit coefficient is established here.'
        ]
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-cycle-size',type=int,default=8)
    parser.add_argument('--max-path-size',type=int,default=8)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=verify(args.max_cycle_size,args.max_path_size)
    text=json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    summary={k:v for k,v in result.items() if k!='window_observables'}
    summary['window_observables']=[{k:v for k,v in row.items()
                                  if k not in ('table_little_endian','nonzero_Walsh_coefficients')}
                                 for row in result['window_observables']]
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
