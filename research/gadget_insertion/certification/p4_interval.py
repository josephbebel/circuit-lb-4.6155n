#!/usr/bin/env python3
"""Rigorous local P4 functional enclosure via interval Taylor integration.

This certifies explicitly supplied covariance points and small neighborhoods,
never the whole distance-type cone.
"""
from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction
import itertools
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'mtree' / 'enclosure'))
from ivarith import IA, asin, PI, dn, up, rsum


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def interval_fraction(value):
    value = Fraction(value)
    near = float(value)
    actual = Fraction.from_float(near)
    return IA(dn(near) if actual > value else near,
              up(near) if actual < value else near)


def exact_endpoints(value, places=15):
    low=Fraction.from_float(float(value.lo));high=Fraction.from_float(float(value.hi))
    scale=10**places
    nlow=(low.numerator*scale)//low.denominator
    nhigh=-((-high.numerator*scale)//high.denominator)
    def decimal(n):
        sign='-' if n<0 else ''
        n=abs(n)
        return f'{sign}{n//scale}.{n%scale:0{places}d}'
    return {'exact_binary64_endpoint_fractions':[str(low),str(high)],
            'outward_decimal_enclosure':[decimal(nlow),decimal(nhigh)]}


class Jet:
    """Bivariate Taylor coefficients through total degree N, in interval lanes."""
    N = 6
    indices = ()
    lookup = {}
    products = ()

    @classmethod
    def configure(cls, degree):
        cls.N = degree
        cls.indices = tuple((a, n-a) for n in range(degree+1) for a in range(n+1))
        cls.lookup = {ab: i for i, ab in enumerate(cls.indices)}
        cls.products = tuple(tuple((cls.lookup[a,b], cls.lookup[x-a,y-b])
                                  for a in range(x+1) for b in range(y+1))
                             for x, y in cls.indices)

    def __init__(self, coefficients):
        self.c = coefficients

    @classmethod
    def constant(cls, value):
        value = IA.of(value)
        return cls([value]+[IA(0.0)]*(len(cls.indices)-1))

    @classmethod
    def variable(cls, value, axis):
        result = cls.constant(value)
        result.c[cls.lookup[(1,0) if axis == 0 else (0,1)]] = IA(1.0)
        return result

    def __add__(self, other):
        if not isinstance(other, Jet):
            return Jet([self.c[0]+other]+self.c[1:])
        return Jet([a+b for a,b in zip(self.c, other.c)])
    __radd__ = __add__

    def __neg__(self):
        return Jet([-x for x in self.c])

    def __sub__(self, other):
        return self + (-other if isinstance(other, Jet) else -IA.of(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return Jet([a*other for a in self.c])
        values = []
        for pairs in self.products:
            total = IA(0.0)
            for i,j in pairs:
                total = total+self.c[i]*other.c[j]
            values.append(total)
        return Jet(values)
    __rmul__ = __mul__

    def reciprocal(self):
        values = [self.c[0].ninv()]
        for index in range(1, len(self.c)):
            total = IA(0.0)
            for i,j in self.products[index]:
                if i:
                    total = total+self.c[i]*values[j]
            values.append(-total*values[0])
        return Jet(values)

    def __truediv__(self, other):
        return self*(other.reciprocal() if isinstance(other, Jet) else IA.of(other).ninv())

    def __rtruediv__(self, other):
        return self.reciprocal()*other

    def sqrt(self):
        values = [self.c[0].nsqrt()]
        denominator = (2*values[0]).ninv()
        for index in range(1, len(self.c)):
            total = self.c[index]
            for i,j in self.products[index]:
                if i and j:
                    total = total-values[i]*values[j]
            values.append(total*denominator)
        return Jet(values)

    def asin(self):
        derivative = (1-self*self).sqrt().reciprocal()
        values = [asin(self.c[0])]
        for a,b in self.indices[1:]:
            axis = 0 if a else 1
            power = a if a else b
            x,y = (a-1,b) if a else (a,b-1)
            total = IA(0.0)
            for i,j in self.products[self.lookup[x,y]]:
                u,v = self.indices[i]
                coefficient = self.c[self.lookup[(u+1,v) if axis == 0 else (u,v+1)]]
                total = total+coefficient*(u+1 if axis == 0 else v+1)*derivative.c[j]
            values.append(total/power)
        return Jet(values)

    def lane_sum(self):
        return Jet([rsum(x, axis=-1) for x in self.c])

    def shift_x_polynomial(self, halfwidth, degree):
        """Exact shift of an x-only polynomial of known degree, no remainder."""
        values = [IA(0.0)]*len(self.c)
        for a in range(degree+1):
            total = IA(0.0)
            for j in range(a,degree+1):
                power = j-a
                bound = halfwidth**power
                delta = IA(-bound,bound) if power%2 else IA(0.0 if power else 1.0,bound)
                total += self.c[self.lookup[j,0]]*math.comb(j,a)*delta
            values[self.lookup[a,0]] = total
        return Jet(values)


Jet.configure(6)
PAIRS = tuple(itertools.combinations(range(6), 2))
SMALL_PAIRS = tuple(itertools.combinations(range(4), 2))
SMALL_LOOKUP = {pair:i for i,pair in enumerate(SMALL_PAIRS)}


def cost_table():
    table = []
    owners = (0,0,1,2,3,3)
    for mask in range(64):
        boundary = [(mask>>i)&1 for i in range(6)]
        values = []
        for internal_mask in range(8):
            incident = [[] for _ in range(4)]
            for i, owner in enumerate(owners):
                incident[owner].append(boundary[i])
            for i in range(3):
                bit = (internal_mask>>i)&1
                incident[i].append(bit); incident[i+1].append(bit)
            values.append(sum(min(row) != max(row) for row in incident))
        table.append(min(values))
    return table


def coefficients():
    table = cost_table()
    result = []
    max_difference = 0
    for i,j in PAIRS:
        rest = tuple(v for v in range(6) if v not in (i,j))
        delta = []
        for mask in range(16):
            base = sum(((mask>>v)&1)<<rest[v] for v in range(4))
            delta.append(table[base|(1<<i)|(1<<j)]-table[base|(1<<i)]
                         -table[base|(1<<j)]+table[base])
        require(max(delta)<=0, 'Boolean submodularity check failed')
        max_difference = max(max_difference, -min(delta))
        numerators = [sum(delta[s]*(-1)**((s&mask).bit_count()) for s in range(16))
                      for mask in range(16)]
        require(all(numerators[m] == 0 for m in range(16) if m.bit_count()%2),
                'Boolean complement-symmetry check failed')
        result.append(numerators)
    return np.asarray(result,dtype=np.int64), Fraction(sum(table),64), max_difference


FOURIER, INDEPENDENT_COST, MAX_DIFFERENCE = coefficients()


def stack_intervals(values):
    return IA(np.asarray([x.lo for x in values]), np.asarray([x.hi for x in values]))


def load_candidate(path, natural_tree=False):
    data = json.loads(path.read_text())
    require(data['path_vertices'] == 4, 'candidate is not P4')
    # Decimal strings are exact rationals; normalization is exact as well.
    types = [tuple(row['type']) for row in data['support']]
    if natural_tree:
        owners = (0,0,1,2,3,3)
        locations = [(i,0) for i in range(4)]+[(i,1) for i in owners]
        tree = {}
        for j,owner in enumerate(owners):
            leaf = 4+j
            distances = [0 if i==leaf else abs(v-owner)+extra+1
                         for i,(v,extra) in enumerate(locations)]
            pattern = tuple(d-distances[0] for d in distances)
            tree[pattern] = Fraction(1,3*2**owner)
        require(set(types)==set(tree),'candidate support does not equal the natural tree face')
        raw = [tree[t] for t in types]
    else:
        raw = [Fraction(str(row['candidate_mass'])) for row in data['support']]
    masses = [w/sum(raw) for w in raw]
    require(all(w>0 for w in masses) and sum(masses)==1, 'invalid rational masses')
    sqrt2 = IA(2.0).sqrt()
    def power(exponent):
        return interval_fraction(Fraction(2)**(exponent//2)) * (sqrt2 if exponent%2 else 1)
    gram = [[sum((interval_fraction(w)*power(-t[i]-t[j]) for t,w in zip(types,masses)),IA(0.0))
             for j in range(10)] for i in range(10)]
    vertex = [[gram[i][j]/(gram[i][i]*gram[j][j]).sqrt() for j in range(10)]
              for i in range(10)]
    owners = (0,0,1,2,3,3)
    edge = [[sum((vertex[x][y] for x in (owners[i],4+i) for y in (owners[j],4+j)),IA(0.0))
             for j in range(6)] for i in range(6)]
    cov = [[IA(1.0) if i==j else edge[i][j]/(edge[i][i]*edge[j][j]).sqrt()
            for j in range(6)] for i in range(6)]
    return cov, [{'type':list(t),'mass':str(w)} for t,w in zip(types,masses)]


def prove_pd(cov, diagonal_shift=0.0):
    """Interval Cholesky proves every matrix in this entry box is PD."""
    n = len(cov)
    L = [[IA(0.0) for _ in range(n)] for _ in range(n)]
    pivots = []
    for i in range(n):
        value = cov[i][i]-diagonal_shift-sum((L[i][k]*L[i][k] for k in range(i)),IA(0.0))
        require(value.lo>0, 'positive definiteness not certified')
        pivots.append([float(value.lo),float(value.hi)])
        L[i][i] = value.sqrt()
        for j in range(i+1,n):
            L[j][i] = (cov[j][i]-sum((L[j][k]*L[i][k] for k in range(i)),IA(0.0)))/L[i][i]
    return pivots


def integrand(cov, rectangles, enclosure, cache=None):
    # Lanes are outer conditioned pairs; first array dimension indexes cells.
    centers = np.asarray([(r[0],r[1]) for r in rectangles])
    radii = np.asarray([(r[2],r[3]) for r in rectangles]) if enclosure else np.zeros_like(centers)
    def tightened(value,label):
        if cache is None:return value
        if not enclosure:
            cache[label]=value
            return value
        center=cache[label]
        # Multivariate mean-value form for each normalized derivative. Work
        # downwards so tighter higher derivatives may tighten lower ones.
        coefficients=list(value.c)
        for index in range(len(coefficients)-1,-1,-1):
            a,b=Jet.indices[index]
            if a+b==Jet.N:continue
            bound=center.c[index]
            for axis,(i,j),factor in [(0,(a+1,b),a+1),(1,(a,b+1),b+1)]:
                delta=IA(-radii[:,axis,None],radii[:,axis,None])
                bound+=factor*coefficients[Jet.lookup[i,j]]*delta
            lower=np.maximum(coefficients[index].lo,bound.lo)
            upper=np.minimum(coefficients[index].hi,bound.hi)
            empty=lower>upper
            coefficients[index]=IA(np.where(empty,np.nan,lower),np.where(empty,np.nan,upper))
        return Jet(coefficients)
    # Form the outer polynomials at the exact midpoint and shift them exactly.
    # Direct repeated interval occurrences of t would destroy cancellations in
    # their small positive principal minors near t=1.
    t = Jet.variable(IA(centers[:,0,None]),0)
    u = Jet.variable(IA((centers[:,1]-radii[:,1])[:,None],
                        (centers[:,1]+radii[:,1])[:,None]),1)
    r = stack_intervals([cov[i][j] for i,j in PAIRS])
    tp = t*r
    t2 = t*t; t3 = t2*t
    om_polynomial = 1-t2*(r*r)
    def shifted(polynomial,degree):
        return polynomial.shift_x_polynomial(radii[:,0,None],degree) if enclosure else polynomial
    om = shifted(om_polynomial,2)
    conditional = []
    for a,b in SMALL_PAIRS:
        params = []
        for i,j in PAIRS:
            rest = [v for v in range(6) if v not in (i,j)]
            k,ell = rest[a],rest[b]
            params.append([cov[k][ell],cov[k][i],cov[k][j],cov[ell][i],cov[ell][j]])
        rkl,rki,rkj,rli,rlj = [stack_intervals([row[z] for row in params]) for z in range(5)]
        numerator = shifted(t*rkl-t2*(rki*rli+rkj*rlj)
                            +t3*(r*(rki*rlj+rkj*rli)-rkl*r*r),3)
        dk = shifted(1-t2*(r*r+rki*rki+rkj*rkj)+t3*(2*r*rki*rkj),3)
        dl = shifted(1-t2*(r*r+rli*rli+rlj*rlj)+t3*(2*r*rli*rlj),3)
        conditional.append(tightened(numerator/(dk*dl).sqrt(),f'conditional_{a}_{b}'))
    expectation = Jet.constant(IA(FOURIER[:,0]/16))
    for pair, value in zip(SMALL_PAIRS, conditional):
        mask = sum(1<<x for x in pair)
        expectation += value.asin()*(2/PI)*IA(FOURIER[:,mask]/16)
    q4 = IA(FOURIER[:,15]/16)
    expectation += conditional[SMALL_LOOKUP[0,1]].asin()*conditional[SMALL_LOOKUP[2,3]].asin()*(4/(PI*PI))*q4
    cross_pairs = ((0,2),(0,3),(1,2),(1,3))
    unscaled = {pair:conditional[index] for index,pair in enumerate(SMALL_PAIRS)}
    def at(i,j):
        return unscaled[tuple(sorted((i,j)))]
    cross_total = Jet.constant(IA(0.0))
    u2=u*u;u3=u2*u
    for i,j in cross_pairs:
        k,ell = [v for v in range(4) if v not in (i,j)]
        r0,z,p,w,v,q = at(i,j),at(k,ell),at(k,i),at(ell,i),at(k,j),at(ell,j)
        denominator = tightened(1-u2*(r0*r0),f'denominator_{i}_{j}')
        numerator = tightened(u*(z-p*w-v*q+r0*p*q)+u3*(r0*v*w-z*r0*r0),f'numerator_{i}_{j}')
        dk = tightened(1-p*p-u2*(r0*r0+v*v-2*r0*p*v),f'dk_{i}_{j}')
        dl = tightened(1-q*q-u2*(r0*r0+w*w-2*r0*q*w),f'dl_{i}_{j}')
        rho = tightened(numerator/(dk*dl).sqrt(),f'rho_{i}_{j}')
        cross_total += conditional[SMALL_LOOKUP[i,j]]*rho.asin()/denominator.sqrt()
    expectation += cross_total*(4/(PI*PI))*q4
    return (expectation*r/(2*PI)/om.sqrt()).lane_sum()


def cell_enclosure(center_coefficients, box_coefficients, lane, hx, hy):
    result = IA(0.0)
    for index,(a,b) in enumerate(Jet.indices):
        degree = a+b
        moment = IA(4.0*hx**(a+1)*hy**(b+1))/((a+1)*(b+1))
        if degree<Jet.N:
            if a%2==0 and b%2==0:
                result += IA(center_coefficients[index].lo[lane],center_coefficients[index].hi[lane])*moment
        elif a%2==0 and b%2==0:
            result += IA(box_coefficients[index].lo[lane],box_coefficients[index].hi[lane])*moment
        else:
            # N is even, so both exponents are odd. The regions where the
            # monomial is positive and negative have equal absolute integrals.
            radius = ((IA(box_coefficients[index].hi[lane])-IA(box_coefficients[index].lo[lane]))*moment/2).hi
            result += IA(-radius,radius)
    return result


def integrate(cov, tolerance, initial=8, batch_size=64, max_cells=100000):
    queue = deque(((i+0.5)/initial,(j+0.5)/initial,1/(2*initial),1/(2*initial))
                  for i in range(initial) for j in range(initial))
    total = IA(0.0)
    accepted = attempted = split = 0
    accepted_area = Fraction(0)
    minimum_radius = 1.0
    started = time.monotonic()
    while queue:
        cells = [queue.popleft() for _ in range(min(batch_size,len(queue)))]
        attempted += len(cells)
        require(attempted<=max_cells, 'cell budget exceeded without certificate')
        with np.errstate(all='ignore'):
            cache={}
            centers = integrand(cov,cells,False,cache).c
            coefficients = integrand(cov,cells,True,cache).c
        for lane, (x,y,hx,hy) in enumerate(cells):
            result = cell_enclosure(centers,coefficients,lane,hx,hy)
            area = 4*hx*hy
            good = bool(np.isfinite(result.lo) and np.isfinite(result.hi)
                        and result.hi-result.lo <= tolerance*area)
            if good:
                total += result
                accepted += 1
                accepted_area += 4*Fraction.from_float(hx)*Fraction.from_float(hy)
                minimum_radius = min(minimum_radius,hx,hy)
            else:
                split += 1
                # Conditioning on an outer pair creates sharper t-dependence
                # than the inner cross-scaling variable. This changes only the
                # subdivision order; both exact children always cover the cell.
                if hx>=hy/8:
                    require(hx>2**-20,'cell became too narrow')
                    queue.extend(((x-hx/2,y,hx/2,hy),(x+hx/2,y,hx/2,hy)))
                else:
                    require(hy>2**-20,'cell became too narrow')
                    queue.extend(((x,y-hy/2,hx,hy/2),(x,y+hy/2,hx,hy/2)))
        print(json.dumps({'attempted':attempted,'accepted':accepted,'remaining':len(queue),
                          'elapsed_seconds':round(time.monotonic()-started,2)}),file=sys.stderr,flush=True)
    require(accepted_area==1,'accepted rectangles did not cover area one exactly')
    return total, {'accepted_cells':accepted,'attempted_cells':attempted,'split_cells':split,
                   'accepted_area_exact':str(accepted_area),
                   'minimum_halfwidth':minimum_radius,'elapsed_seconds':time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,default=Path(__file__).parents[1]/'cone_search_P4.json')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--tree',action='store_true',help='certify exact natural tree masses instead')
    parser.add_argument('--degree',type=int,default=6)
    parser.add_argument('--tolerance',type=float,default=2e-7)
    parser.add_argument('--initial',type=int,default=8)
    parser.add_argument('--batch-size',type=int,default=64)
    parser.add_argument('--max-cells',type=int,default=100000)
    args = parser.parse_args()
    require(4<=args.degree<=12 and args.degree%2==0,'Taylor degree must be even, between four and twelve')
    require(0<args.initial<=1024 and args.initial&(args.initial-1)==0,
            'initial grid must be a power of two at most 1024')
    require(math.isfinite(args.tolerance) and args.tolerance>0,'integration tolerance must be positive and finite')
    require(args.batch_size>0 and args.max_cells>0,'positive batch size and cell budget required')
    from mpmath import iv
    require(bool(iv.mpf(float(PI.lo))<iv.pi) and bool(iv.pi<iv.mpf(float(PI.hi))),
            'pi interval did not pass independent high-precision check')
    Jet.configure(args.degree)
    cov,masses = load_candidate(args.candidate,args.tree)
    sys.path.insert(0,str(Path(__file__).parents[1]))
    from path_cone import enumerate_types
    type_set=set(enumerate_types(4))
    require(len(type_set)==4655 and all(len(t)==10 and max(map(abs,t))<=4 for t in type_set),
            'P4 cone type or generator-size check failed')
    require(all(tuple(row['type']) in type_set for row in masses),'candidate is outside the declared cone')
    arithmetic_source=(ROOT/'mtree'/'enclosure'/'ivarith.py').read_bytes()
    arithmetic_blob=hashlib.sha1(b'blob '+str(len(arithmetic_source)).encode()+b'\0'+arithmetic_source).hexdigest()
    require(arithmetic_blob=='b99d01f757f2b3d75151ba9403a547c6b8c60477',
            'interval arithmetic source differs from the audited upstream version')
    pivots = prove_pd(cov)
    integral,stats = integrate(cov,args.tolerance,args.initial,args.batch_size,args.max_cells)
    value = (interval_fraction(INDEPENDENT_COST)+integral)/4
    radius = Fraction(1,2**20)
    radius_iv = interval_fraction(radius)
    local_pivots = prove_pd(cov,float(5*radius))
    corr_limit = interval_fraction(Fraction(24,25))
    require(all((cov[i][j]+radius_iv).hi<corr_limit.lo
                and (-cov[i][j]+radius_iv).hi<corr_limit.lo for i,j in PAIRS),
            'local correlation box exceeds the derivative-bound range')
    require(MAX_DIFFERENCE<=2 and PI.lo>3,'local Lipschitz prerequisites failed')
    # |d Phi / d r_ij| <= 1/(pi sqrt(1-(24/25)^2)) < 25/21.
    # Fifteen correlation coordinates and division by four vertices give
    # an upper per-vertex change of (125/28)*radius.
    local_upper = value+interval_fraction(Fraction(125,28))*radius_iv
    require(Fraction.from_float(float(local_upper.hi))<Fraction(6807,50000),
            'the requested local upper bound 0.13614 was not certified')
    exact_gram_diagonal=[sum((Fraction(row['mass'])*Fraction(2)**(-row['type'][i])
                             for row in masses),Fraction(0)) for i in range(10)]
    require(min(exact_gram_diagonal)>=Fraction(4,5),'mass-neighborhood diagonal bound failed')
    mass_radius=Fraction(1,2**30)
    require(640*mass_radius<=radius and 160*mass_radius<=Fraction(1,2),
            'mass-neighborhood normalization estimate failed')
    record = {'status':'rigorous local point enclosure; not a whole-cone bound',
              'point':'natural_tree' if args.tree else 'rationalized_nonuniform_candidate',
              'rational_masses':masses,'covariance_lo':[[float(x.lo) for x in row]for row in cov],
              'covariance_hi':[[float(x.hi) for x in row]for row in cov],
              'cholesky_pivots':pivots,'per_vertex_enclosure':[float(value.lo),float(value.hi)],
              'per_vertex_exact_endpoints':exact_endpoints(value),
              'local_correlation_box_radius':str(radius),
              'local_box_shifted_cholesky_pivots':local_pivots,
              'local_box_per_vertex_upper_bound':float(local_upper.hi),
              'local_box_upper_endpoint_fraction':str(Fraction.from_float(float(local_upper.hi))),
              'local_box_rational_upper_bound':'6807/50000',
              'local_simplex_L1_radius':str(mass_radius),
              'exact_gram_diagonal':[str(x)for x in exact_gram_diagonal],
              'scope':'Point, correlation box, and stated L1 simplex neighborhood only; not the entire cone.',
              'upstream_interval_arithmetic_blob':arithmetic_blob,
              'integration':stats,'taylor_degree':args.degree,'max_boolean_mixed_difference':MAX_DIFFERENCE}
    text=json.dumps(record,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()
