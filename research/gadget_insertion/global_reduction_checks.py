#!/usr/bin/env python3
"""Exact finite checks supporting GLOBAL_REDUCTION.md.

These exercise graph identities and layout inequalities. The general theorem
is proved in the note; no finite list of graphs replaces that proof.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import itertools
import json
from pathlib import Path
import random

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


def require(test, message):
    if not test:
        raise ArithmeticError(message)


@dataclass(frozen=True)
class Graph:
    n: int
    edges: tuple[tuple[int, int], ...]

    def adjacency(self):
        rows = [[] for _ in range(self.n)]
        for e,(u,v) in enumerate(self.edges):
            rows[u].append((v,e)); rows[v].append((u,e))
        return rows

    def component(self, start, skip=()):
        skip=set(skip); rows=self.adjacency()
        found={start}; stack=[start]
        while stack:
            u=stack.pop()
            for v,e in rows[u]:
                if e not in skip and v not in found:
                    found.add(v);stack.append(v)
        return frozenset(found)

    def bridges(self, omitted=None):
        rows=self.adjacency(); depth=[-1]*self.n; low=[0]*self.n
        result=[]
        def visit(u,parent_edge,d):
            depth[u]=low[u]=d
            for v,e in rows[u]:
                if e==omitted or e==parent_edge:
                    continue
                if depth[v]<0:
                    visit(v,e,d+1)
                    low[u]=min(low[u],low[v])
                    if low[v]>depth[u]:result.append(e)
                else:low[u]=min(low[u],depth[v])
        for u in range(self.n):
            if depth[u]<0:visit(u,None,0)
        return result

    def check_class(self, degree=3):
        require(self.n>=2,'too few vertices')
        require(all(0<=u<self.n and 0<=v<self.n and u!=v for u,v in self.edges),
                'not loopless')
        require(all(len(row)==degree for row in self.adjacency()),'wrong degree')
        require(len(self.component(0))==self.n,'disconnected')
        require(not self.bridges(),'bridge found')

    def two_cuts(self):
        result={}; seen=set(); all_vertices=frozenset(range(self.n))
        for e in range(len(self.edges)):
            for f in self.bridges(omitted=e):
                cut=tuple(sorted((e,f)))
                if cut in seen:continue
                seen.add(cut)
                shore=self.component(self.edges[f][0],cut)
                other=all_vertices-shore
                require(shore and other,'reported 2-cut does not disconnect')
                canonical=min((shore,other),key=lambda s:(len(s),tuple(sorted(s))))
                result[canonical]=cut
        return result

    def width(self,order):
        require(len(order)==self.n and set(order)==set(range(self.n)),'invalid layout')
        rows=self.adjacency(); processed=set(); crossing=maximum=0
        for u in order:
            crossing+=len(rows[u])-2*sum(v in processed for v,_ in rows[u])
            processed.add(u);maximum=max(maximum,crossing)
        require(crossing==0,'layout edge accounting failed')
        return maximum

    def simple_edges(self):
        return sorted({tuple(sorted(e)) for e in self.edges})

    def bfs_order(self):
        rows=self.adjacency();order=[0];seen={0}
        for u in order:
            for v,_ in rows[u]:
                if v not in seen:seen.add(v);order.append(v)
        require(len(order)==self.n,'BFS input disconnected')
        return order


def matching_witness(g):
    """A solver proposes edges; all used matching claims are checked exactly.

    Optimality certificates from the floating-point MILP solver are not
    trusted. Callers verify the required cardinality, whose upper bound is
    independently known from vertex count or an explicit Tutte barrier.
    """
    edges=g.simple_edges();a=np.zeros((g.n,len(edges)))
    for e,(u,v) in enumerate(edges):a[u,e]=a[v,e]=1
    result=milp(-np.ones(len(edges)),integrality=np.ones(len(edges)),
                bounds=Bounds(0,1),constraints=LinearConstraint(a,0,1),
                options={'time_limit':30})
    require(result.x is not None,'no matching proposal')
    chosen=[edge for edge,x in zip(edges,result.x) if x>0.5]
    require(len({v for edge in chosen for v in edge})==2*len(chosen),
            'proposed matching is not vertex-disjoint')
    require(all(edge in edges for edge in chosen),'proposed nonedge')
    return chosen


def generalized_petersen(n,k):
    return Graph(2*n,tuple([(i,(i+1)%n) for i in range(n)]+
                           [(i,n+i) for i in range(n)]+
                           [(n+i,n+(i+k)%n) for i in range(n)]))


def random_cubic(n,seed):
    rng=random.Random(seed)
    for _ in range(10000):
        stubs=[v for v in range(n) for _ in range(3)];rng.shuffle(stubs)
        edges=[tuple(sorted(stubs[i:i+2])) for i in range(0,3*n,2)]
        if any(u==v for u,v in edges) or len(set(edges))<len(edges):continue
        g=Graph(n,tuple(edges))
        if len(g.component(0))==n and not g.bridges():return g
    raise ArithmeticError('random family construction failed')


def completed_shore(g,shore):
    vertices=sorted(shore);loc={v:i for i,v in enumerate(vertices)}
    boundary=[];edges=[]
    for u,v in g.edges:
        if u in shore and v in shore:edges.append((loc[u],loc[v]))
        elif u in shore:boundary.append(u)
        elif v in shore:boundary.append(v)
    require(len(boundary)==2 and boundary[0]!=boundary[1],'invalid shore boundary')
    edges.append(tuple(loc[u] for u in boundary))
    answer=Graph(len(vertices),tuple(edges));answer.check_class()
    return answer,vertices


def factor_on_three_edge_connected(g):
    g.check_class();require(not g.two_cuts(),'not 3-edge-connected')
    if g.n==2:return {'order':2,'bounded_triple_edge_exception':True}
    simple=set(g.simple_edges())
    require(len(simple)==len(g.edges),'nontrivial core is not simple')
    matching=matching_witness(g)
    require(2*len(matching)==g.n,'Petersen matching missing')
    owner={v:i for i,edge in enumerate(matching) for v in edge}
    matchset={frozenset(e) for e in matching};qe=[];representative={}
    for u,v in g.edges:
        if frozenset((u,v)) in matchset:continue
        a,b=owner[u],owner[v]
        require(a!=b,'loop after contraction')
        qe.append((a,b));representative.setdefault(frozenset((a,b)),(u,v))
    q=Graph(len(matching),tuple(qe));q.check_class(4)
    require(not q.two_cuts(),'contracted core has a 2-edge cut')
    packing=matching_witness(q)
    require(len(packing)==q.n//2,'factor leaves more than one matching pair')
    paths=[];used=set()
    for a,b in packing:
        u,v=representative[frozenset((a,b))]
        uu=next(x for x in matching[owner[u]] if x!=u)
        vv=next(x for x in matching[owner[v]] if x!=v)
        path=(uu,u,v,vv)
        require(len(set(path))==4 and not (set(path)&used),'overlapping P4')
        require(all(tuple(sorted((x,y))) in simple for x,y in zip(path,path[1:])),'invalid P4')
        paths.append(path);used.update(path)
    leftover=set(range(g.n))-used
    require(len(leftover) in (0,2),'wrong residual size')
    if leftover:require(tuple(sorted(leftover)) in simple,'residual is not P2')
    return {'order':g.n,'P4_paths':len(paths),'P2_paths':int(bool(leftover)),
            'all_vertices_covered':True}


def verify_small_shore_core(g,T,cuts):
    require(g.n>3*T,'small graph outside lemma regime')
    require(all(min(len(a),g.n-len(a))<T for a in cuts),'balanced cut present')
    all_vertices=frozenset(range(g.n));small=set()
    for a in cuts:
        for shore in (a,all_vertices-a):
            if len(shore)<T:small.add(shore)
    maximal=sorted((a for a in small if not any(a<b for b in small)),
                   key=lambda a:tuple(sorted(a)))
    for a,b in itertools.combinations(maximal,2):
        require(not a&b,'maximal shores overlap')
        require(not any((u in a and v in b) or (u in b and v in a) for u,v in g.edges),
                'maximal shores adjacent')
    removed=set().union(*maximal) if maximal else set()
    vertices=sorted(set(range(g.n))-removed);loc={v:i for i,v in enumerate(vertices)}
    edges=[(loc[u],loc[v]) for u,v in g.edges if u in loc and v in loc]
    expansions=[]
    for a in maximal:
        outside=[]
        for u,v in g.edges:
            if u in a and v not in a:outside.append(v)
            elif v in a and u not in a:outside.append(u)
        require(len(outside)==2 and len(set(outside))==2,'bad external endpoints')
        require(all(v in loc for v in outside),'attachment outside core')
        edges.append(tuple(loc[v] for v in outside));expansions.append((a,tuple(outside)))
    core=Graph(len(vertices),tuple(edges));core.check_class()
    require(not core.two_cuts(),'replacement core still has a 2-cut')
    orders=[list(range(core.n)),list(reversed(range(core.n)))]
    if core.n>2:orders.append(core.bfs_order())
    checked=[]
    for order in orders:
        original_order=[vertices[i] for i in order];position={v:i for i,v in enumerate(original_order)}
        after={v:[] for v in original_order}
        for a,(x,y) in expansions:
            earlier=min((x,y),key=position.__getitem__)
            after[earlier].append(sorted(a))
        expanded=[]
        for v in original_order:
            expanded.append(v)
            for a in after[v]:expanded.extend(a)
        w0=core.width(order);w1=g.width(expanded)
        require(2*(w1-w0)<=3*T,'expansion inequality failed')
        checked.append({'core_width':w0,'expanded_width':w1})
    return {'T':T,'maximal_small_shores':len(maximal),'core_order':core.n,
            'largest_shore':max(map(len,maximal),default=0),'layouts':checked}


def two_sum(g,h,e=0,f=0):
    a,b=g.edges[e];c,d=h.edges[f]
    edges=[edge for i,edge in enumerate(g.edges) if i!=e]
    edges.extend((u+g.n,v+g.n) for i,(u,v) in enumerate(h.edges) if i!=f)
    edges.extend(((a,c+g.n),(b,d+g.n)))
    answer=Graph(g.n+h.n,tuple(edges));answer.check_class()
    return answer


def sharp_lift(s):
    # A doubled cycle is represented by its two opposite directed cycles.
    directed=[];n=s
    for u in range(s):
        for v in ((u+1)%s,(u-1)%s):
            a,b,c=n,n+1,n+2;n+=3
            directed.extend(((a,b),(b,a),(a,c),(c,a),(b,c),(u,b),(c,v)))
    q=Graph(n,tuple(directed));q.check_class(4)
    matching=matching_witness(q)
    require(n==7*s and len(matching)==3*s,'sharp matching count failed')
    residual=set(range(s,n));sizes=[];rows=q.adjacency()
    while residual:
        start=next(iter(residual));found={start};stack=[start]
        while stack:
            u=stack.pop()
            for v,_ in rows[u]:
                if v>=s and v not in found:found.add(v);stack.append(v)
        residual-=found;sizes.append(len(found))
    sizes.sort()
    require(sizes==[3]*(2*s),'wrong Tutte barrier components')
    h=Graph(2*n,tuple([(2*v,2*v+1) for v in range(n)]+
                      [(2*u,2*v+1) for u,v in directed]))
    h.check_class()
    require(len({frozenset(e) for e in h.edges})==len(h.edges),'lift not simple')
    require(all((u-v)%2 for u,v in h.edges),'lift not bipartite')
    return h,{'base_order':s,'contracted_order':n,'maximum_matching_size':len(matching),
              'Tutte_barrier_size':s,'odd_components_after_deletion':2*s,
              'simple_bipartite_cubic_lift_order':h.n}


def verify():
    families={
        'triple_edge':Graph(2,((0,1),)*3),
        'K4':Graph(4,tuple(itertools.combinations(range(4),2))),
        'K33':Graph(6,tuple((u,v) for u in range(3) for v in range(3,6))),
        'cube':Graph(8,tuple((u,u^(1<<b)) for u in range(8) for b in range(3) if u<u^(1<<b))),
        'Petersen':generalized_petersen(5,2),
        'Heawood':Graph(14,tuple((j,7+(j+d)%7) for j in range(7) for d in (0,1,3))),
        'dodecahedral':generalized_petersen(10,2),
    }
    for k in range(3,13):
        families[f'prism_{k}']=generalized_petersen(k,1)
    for n in (8,10,14,20,30,40):
        for seed in (18,37,91):
            g=random_cubic(n,seed)
            if len(g.component(0))==n and not g.bridges():families[f'random_{n}_{seed}']=g
    atom=families['cube'];chain=atom
    for index in range(1,9):
        chain=two_sum(chain,atom,e=len(chain.edges)-3)
        if index in (1,3,8):families[f'cube_chain_{index+1}']=chain
    families['parallel_shore']=two_sum(families['cube'],families['triple_edge'])
    core=families['cube'];expanded=core
    for edge in core.edges:
        e=next(i for i,current in enumerate(expanded.edges) if frozenset(current)==frozenset(edge))
        expanded=two_sum(expanded,families['K33'],e=e)
    families['cube_all_edges_expanded']=expanded
    sharp=[]
    for s in (3,4,7):
        g,record=sharp_lift(s);families[f'sharp_lift_{s}']=g;sharp.append(record)

    cases=[];total_cuts=0;core_checks=0;factor_checks=0
    for name,g in families.items():
        g.check_class();cuts=g.two_cuts();total_cuts+=len(cuts)
        row={'name':name,'order':g.n,'two_edge_cuts':len(cuts)}
        if not cuts:
            row['factor']=factor_on_three_edge_connected(g);factor_checks+=1
        # Every cut completion and a concrete concatenation are checked.
        for a in cuts:
            b=frozenset(range(g.n))-a
            ga,va=completed_shore(g,a);gb,vb=completed_shore(g,b)
            combined=va+vb
            require(g.width(combined)<=max(ga.width(range(ga.n)),gb.width(range(gb.n)))+2,
                    'concatenation inequality failed')
        local=[]
        for T in (3,5,8,15):
            if g.n>3*T and all(min(len(a),g.n-len(a))<T for a in cuts):
                local.append(verify_small_shore_core(g,T,cuts));core_checks+=1
        if local:row['small_shore_core_checks']=local
        cases.append(row)
    require(core_checks>0 and factor_checks>0 and total_cuts>0,'missing coverage')
    return {'status':'EXACT GRAPH AND LAYOUT CHECKS PASSED; SEE GENERAL PROOF',
            'graph_cases':len(cases),'two_cut_completions_checked':2*total_cuts,
            'three_edge_connected_factor_cases':factor_checks,
            'small_shore_core_parameter_cases':core_checks,
            'sharp_fixed_matching_examples':sharp,'cases':cases,
            'scope':['Finite implementation checks, not a substitute for the general proof.',
                     'MILP proposes matchings; vertex-disjointness and required cardinality are checked on exact integer graphs.',
                     'Floating-point solver optimality reports are not used as proof inputs.',
                     'No Gaussian cone maximum or leaderboard target is certified by this script.']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();data=verify();text=json.dumps(data,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':main()
