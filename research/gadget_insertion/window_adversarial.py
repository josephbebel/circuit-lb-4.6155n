#!/usr/bin/env python3
"""Limited numerical H7 searches on the projected P8 tree face.

This is a feasibility/adversarial probe, not a global search or certificate.
The two objectives are forward H7 and its average with the reversed window.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr,ndtri
from scipy.stats import qmc

from path_cone import covariance,generators,tree_support


def tables():
    forward=[]
    for mask in range(256):
        b=[(mask>>i)&1 for i in range(8)]
        d=[b[i]!=b[i+1] for i in range(7)]
        chain=1;value=0
        for j in range(7):
            chain*=int(d[6-j]);value+=(-1)**j*chain
        forward.append(value)
    forward=np.asarray(forward,dtype=float)
    reversal=np.asarray([sum(((mask>>j)&1)<<(7-j) for j in range(8)) for mask in range(256)])
    reverse=forward[reversal]
    if not np.array_equal(forward,forward[255-np.arange(256)]):
        raise ArithmeticError('complement symmetry failed')
    return np.column_stack((forward,(forward+reverse)/2))


def normal_points(power,seed):
    u=qmc.Sobol(8,scramble=True,seed=seed).random_base2(power)
    return ndtri(np.clip(u,np.finfo(float).eps,1-np.finfo(float).eps))


def conditional_samples(cov,table,normals):
    ones=np.ones(8)
    alpha=.95/float(ones@np.linalg.solve(cov,ones))
    residual=cov-alpha*np.ones_like(cov)
    factor=np.linalg.cholesky((residual+residual.T)/2)
    thresholds=-(normals@factor.T)/math.sqrt(alpha)
    order=np.argsort(thresholds,axis=1,kind='stable')
    sorted_thresholds=np.take_along_axis(thresholds,order,axis=1)
    masks=255^np.bitwise_xor.accumulate(np.left_shift(np.uint64(1),order.astype(np.uint64)),axis=1)
    weights=np.diff(ndtr(sorted_thresholds),axis=1)
    # Below/above all crossings every bit agrees and both objectives vanish.
    return np.sum(table[masks[:,:-1]]*weights[:,:,None],axis=1)


def run(train_power=13,holdout_power=18,repeats=6):
    started=time.monotonic()
    support,exact_natural=tree_support(8)
    natural=np.asarray(exact_natural,dtype=float)
    gen=generators(support)
    chosen=np.arange(1,9)
    table=tables()
    training=normal_points(train_power,99117)
    def masses(z):
        w=natural*np.exp(np.r_[z,0]);return w/w.sum()
    def cov_of(w):
        full=covariance(8,w,gen);return full[np.ix_(chosen,chosen)]
    def estimate(z):
        return conditional_samples(cov_of(masses(z)),table,training).mean(axis=0)
    rng=np.random.default_rng(99119)
    probes=[np.zeros(9)]
    probes.extend(np.clip(rng.normal(0,.35,size=(64,9)),-.8,.8))
    for j in range(9):
        for sign in (-1,1):
            z=np.zeros(9);z[j]=sign*.35;probes.append(z)
    evaluations=np.asarray([estimate(z) for z in probes])
    candidates=[{'name':'natural','z':np.zeros(9),'training_values':evaluations[0]}]
    fits=[]
    for objective in range(2):
        best=int(np.argmax(evaluations[:,objective]))
        for start_name,z0 in [('natural',np.zeros(9)),('best_probe',probes[best])]:
            result=minimize(lambda z:-estimate(z)[objective],z0,method='L-BFGS-B',
                            bounds=[(-.8,.8)]*9,
                            options={'maxiter':80,'ftol':1e-11,'gtol':1e-6,'eps':1e-5})
            vals=estimate(result.x)
            name=('forward' if objective==0 else 'symmetric')+'_'+start_name
            fits.append({'name':name,'success':bool(result.success),'message':str(result.message),
                         'iterations':int(result.nit),'function_evaluations':int(result.nfev),
                         'training_values':vals.tolist(),'log_mass_coordinates':result.x.tolist()})
            candidates.append({'name':name,'z':result.x,'training_values':vals})
    matrices=[cov_of(masses(c['z'])) for c in candidates]
    holdouts=[[] for _ in candidates]
    for repeat in range(repeats):
        normals=normal_points(holdout_power,99231+repeat)
        for i,cov in enumerate(matrices):
            holdouts[i].append(conditional_samples(cov,table,normals).mean(axis=0).tolist())
        print(json.dumps({'holdout':repeat+1,'seconds':round(time.monotonic()-started,2)}),flush=True)
    natural_holdout=np.asarray(holdouts[0])
    rows=[]
    for candidate,cov,values in zip(candidates,matrices,holdouts):
        values=np.asarray(values);paired=values-natural_holdout
        w=masses(candidate['z'])
        rows.append({
            'name':candidate['name'],'masses':w.tolist(),
            'L1_distance_from_natural':float(np.abs(w-natural).sum()),
            'minimum_covariance_eigenvalue':float(np.linalg.eigvalsh(cov)[0]),
            'training_values_forward_symmetric':candidate['training_values'].tolist(),
            'holdout_mean_forward_symmetric':values.mean(axis=0).tolist(),
            'holdout_standard_error_heuristic':(values.std(axis=0,ddof=1)/np.sqrt(repeats)).tolist(),
            'paired_gain_over_natural':paired.mean(axis=0).tolist(),
            'paired_gain_standard_error_heuristic':(paired.std(axis=0,ddof=1)/np.sqrt(repeats)).tolist(),
            'holdout_values':values.tolist(),
            'projected_covariance':cov.tolist(),
        })
    return {
        'status':'LIMITED NUMERICAL ADVERSARIAL FEASIBILITY SEARCH; NO GLOBAL OR INTERVAL CERTIFICATE',
        'objective_order':['forward_H7','half_forward_plus_reversed_H7'],
        'path_core_vertices':8,'tree_face_types':[list(t)for t in support],
        'natural_masses_exact':[str(w)for w in exact_natural],
        'projection_boundary_indices':chosen.tolist(),
        'projection_owners':list(range(8)),
        'relative_mass_log_coordinate_bounds':[-.8,.8],
        'last_log_coordinate_fixed':0,
        'probe_count':len(probes),'training_power':train_power,'training_seed':99117,
        'random_probe_seed':99119,'holdout_power':holdout_power,'holdout_repeats':repeats,
        'holdout_seeds':list(range(99231,99231+repeats)),
        'method':'positive common-factor conditional integration with a smooth Cholesky residual; arbitrary exact 256-entry Boolean tables',
        'optimizer_results':fits,'results':rows,'elapsed_seconds':time.monotonic()-started,
        'scope':'Only the ten-type tree face and bounded relative log-mass perturbations. Any above-target candidate is numerical evidence until its integration error is certified. No uniform upper bound is established.'
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=run()
    text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':main()
