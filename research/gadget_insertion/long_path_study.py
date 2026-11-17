#!/usr/bin/env python3
"""Direct Gaussian/linear-DP natural-tree study for long path gadgets.

This uses no exponential Boolean cost table. Common Gaussian tree samples
couple all prefix lengths. The exact-mean sum of old P2 costs is a control
variate; randomized Sobol errors remain heuristic, not certificates.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.special import ndtri
from scipy.stats import qmc

from path_gadgets import common_tree_covariance, cost, prefix_indices, tree_covariance


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def study(sizes, power, repeats, seed, batch_power):
    sizes=sorted(set(sizes))
    require(all(k>=2 and k%2==0 for k in sizes),'sizes must be positive even path lengths')
    max_size=max(sizes)
    covariance=common_tree_covariance(max_size)
    eigenvalues,eigenvectors=np.linalg.eigh(covariance)
    require(eigenvalues[0]>-1e-10,'common numerical covariance is not PSD')
    keep=eigenvalues>1e-11
    require(int(keep.sum())==max_size+2,'unexpected tree-wave rank')
    factor=eigenvectors[:,keep][:,::-1]*np.sqrt(eigenvalues[keep][::-1])[None,:]
    reconstructed=factor@factor.T
    reconstruction_error=float(np.max(np.abs(reconstructed-covariance)))
    require(reconstruction_error<1e-10,'factorization reconstruction error too large')
    indices={k:prefix_indices(k) for k in sizes}
    marginal_errors={k:float(np.max(np.abs(covariance[np.ix_(idx,idx)]-tree_covariance(k))))
                     for k,idx in indices.items()}
    require(max(marginal_errors.values())<1e-12,'prefix covariance mismatch')

    tree_path=Path(__file__).parents[1]/'mtree_refinement'/'certificate.json'
    record=json.loads(tree_path.read_text())['M_tree']
    tree_lower=Fraction(int(record['lower_numerator']),int(record['denominator']))
    tree_upper=Fraction(int(record['upper_numerator']),int(record['denominator']))
    baseline=float((tree_lower+tree_upper)/4)
    baseline_exact_enclosure=[str(tree_lower/2),str(tree_upper/2)]
    use_batch=min(batch_power,power)
    samples=2**power
    estimates={k:[] for k in sizes}
    raw_estimates={k:[] for k in sizes}
    gain_estimates={k:[] for k in sizes}
    windows=[p for p in (3,5,7,9) if p<max_size]
    window_estimates={p:[] for p in windows}
    window_change_control={p:[] for p in windows}
    window_raw={p:[] for p in windows}
    change_mean=math.acos((2+math.sqrt(2))/4)/math.pi
    started=time.monotonic()
    for repeat in range(repeats):
        sampler=qmc.Sobol(factor.shape[1],scramble=True,seed=seed+repeat)
        sums={k:0 for k in sizes}
        gain_sums={k:0 for k in sizes}
        window_sums={p:0 for p in windows}
        window_D_sums={p:0 for p in windows}
        pair_total_sum=0
        for _ in range(2**(power-use_batch)):
            uniform=sampler.random(2**use_batch)
            normals=ndtri(np.clip(uniform,np.finfo(float).eps,1-np.finfo(float).eps))
            bits=(normals@factor.T)<0
            # Pair j contains vertices 2j and 2j+1. Its four fixed boundary
            # scores are the entering backbone edge, two pendant edges,
            # and the exiting backbone edge. The pair's own backbone edge
            # is optimized. F2 is evaluated by its exact Boolean formula.
            a=bits[:,0:2*max_size:4]
            b=bits[:,1:2*max_size:4]
            c=bits[:,3:2*max_size:4]
            d=bits[:,4:2*max_size+1:4]
            pair_costs=(a!=b).astype(np.int32)+(c!=d).astype(np.int32)
            pair_costs+=((a==b)&(c==d)&(a!=c)).astype(np.int32)
            pair_prefix=np.cumsum(pair_costs,axis=1,dtype=np.int32)
            pair_total_sum+=int(np.sum(pair_prefix[:,-1],dtype=np.int64))
            pendant=bits[:,1:2*max_size:2]
            changes=pendant[:,1:]!=pendant[:,:-1]
            for p in windows:
                chain=changes[:,p-1:].copy()
                h=chain.astype(np.int32)
                for j in range(1,p):
                    chain&=changes[:,p-1-j:max_size-1-j]
                    h+=((-1)**j)*chain.astype(np.int32)
                require(bool(np.all((h==0)|(h==1))),'odd-window Boolean range failed')
                window_sums[p]+=int(np.sum(h,dtype=np.int64))
                window_D_sums[p]+=int(np.sum(changes[:,p-1:],dtype=np.int64))
            for k,idx in indices.items():
                values=cost(bits[:,idx])
                separate=pair_prefix[:,k//2-1]
                gain=separate-values
                require(bool(np.all(gain>=0)),'joint optimization lost against separate pairs')
                sums[k]+=int(np.sum(values,dtype=np.int64))
                gain_sums[k]+=int(np.sum(gain,dtype=np.int64))
        for k in sizes:
            raw_estimates[k].append(sums[k]/(samples*k))
            gain_estimates[k].append(gain_sums[k]/(samples*k))
            estimates[k].append(baseline-gain_estimates[k][-1])
        for p in windows:
            raw=window_sums[p]/(samples*(max_size-p))
            window_raw[p].append(raw)
            window_estimates[p].append(baseline+raw-pair_total_sum/(samples*max_size))
            window_change_control[p].append(change_mean+
                (window_sums[p]-window_D_sums[p])/(samples*(max_size-p)))
        print(json.dumps({'completed_scramble':repeat+1,'scrambles':repeats,
                          'seconds':round(time.monotonic()-started,2),
                          'control_variate_estimates':{str(k):estimates[k][-1] for k in sizes},
                          'odd_window_P2_control_estimates':{str(p):window_estimates[p][-1] for p in windows}}),flush=True)
    rows=[]
    target=Fraction(1341,10000)
    for k in sizes:
        values=np.asarray(estimates[k])
        raw=np.asarray(raw_estimates[k])
        mean=float(values.mean())
        rows.append({
            'path_vertices':k,
            'mean_tree_cost_per_vertex':mean,
            'replicate_control_variate_values':values.tolist(),
            'standard_error_across_scrambles_heuristic':float(values.std(ddof=1)/np.sqrt(repeats)),
            'raw_direct_mean':float(raw.mean()),
            'raw_direct_replicate_values':raw.tolist(),
            'raw_direct_standard_error_heuristic':float(raw.std(ddof=1)/np.sqrt(repeats)),
            'mean_gain_over_separate_P2_per_vertex':float(np.mean(gain_estimates[k])),
            'difference_from_target_0_1341':mean-float(target),
            'prefix_covariance_max_absolute_error':marginal_errors[k],
            'circuit_coefficient_if_this_tree_value_were_universal':1+1/(2*mean),
        })
    window_rows=[]
    for p in windows:
        vals=np.asarray(window_estimates[p])
        raw=np.asarray(window_raw[p])
        dc=np.asarray(window_change_control[p])
        window_rows.append({
            'odd_window_length':p,'matching_score_coordinates':p+1,
            'positions_per_sample':max_size-p,
            'mean_with_P2_control':float(vals.mean()),
            'replicate_P2_control_values':vals.tolist(),
            'P2_control_standard_error_heuristic':float(vals.std(ddof=1)/np.sqrt(repeats)),
            'mean_with_sign_change_control':float(dc.mean()),
            'sign_change_control_standard_error_heuristic':float(dc.std(ddof=1)/np.sqrt(repeats)),
            'raw_mean':float(raw.mean()),
            'raw_standard_error_heuristic':float(raw.std(ddof=1)/np.sqrt(repeats)),
            'difference_from_target_0_1341':float(vals.mean())-float(target),
        })
    return {
        'status':'EXPLORATORY NATURAL-TREE NUMERICS, NOT A GLOBAL OR INTERVAL CERTIFICATE',
        'method':'direct common-tree Gaussian scrambled Sobol sampling; exact O(k) Boolean dynamic programming; exact-mean P2 control variate',
        'sizes':sizes,'points_per_scramble':samples,'scrambles':repeats,
        'first_seed':seed,'seeds':list(range(seed,seed+repeats)),
        'batch_power':use_batch,'common_covariance_dimension':len(covariance),
        'gaussian_rank':factor.shape[1],
        'factorization_max_absolute_error':reconstruction_error,
        'minimum_eigenvalue_binary64_diagnostic':float(eigenvalues[0]),
        'positive_eigenvalue_threshold':1e-11,
        'P2_per_vertex_control_mean_exact_enclosure':baseline_exact_enclosure,
        'reported_frontier_cutwidth_target':'1341/10000',
        'uncertainty':'Replicate variation is heuristic; it is not a rigorous error bound or confidence guarantee.',
        'scope':'Natural-tree points only. Their values constrain what a universal path-gadget upper bound could achieve; no global maximum, limiting theorem, or record is established.',
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'elapsed_seconds':time.monotonic()-started,'results':rows,
        'odd_window_definition':'H_p(i)=sum_{j=0}^{p-1}(-1)^j product_{a=0}^j D_(i-a), with D_i=1[pendant_sign_i differs from pendant_sign_(i+1)]',
        'odd_window_sign_change_control_mean':change_mean,
        'odd_window_results':window_rows,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes',type=int,nargs='+',default=[16,24,32,48,64,128])
    parser.add_argument('--power',type=int,default=20)
    parser.add_argument('--repeats',type=int,default=8)
    parser.add_argument('--seed',type=int,default=88123)
    parser.add_argument('--batch-power',type=int,default=13)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    require(4<=args.power<=24 and 2<=args.repeats<=64,'invalid study size')
    require(4<=args.batch_power<=18 and max(args.sizes)<=512,'invalid batch/path size')
    result=study(args.sizes,args.power,args.repeats,args.seed,args.batch_power)
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()
