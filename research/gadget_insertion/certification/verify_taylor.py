#!/usr/bin/env python3
"""Independent exact-rational checks of the interval Taylor engine."""
from fractions import Fraction as F
import json
import math
import numpy as np

from p4_interval import IA, Jet, cell_enclosure, require


def enclosed(interval, exact):
    low = F.from_float(float(np.asarray(interval.lo).item()))
    high = F.from_float(float(np.asarray(interval.hi).item()))
    require(low<=exact<=high, f'exact coefficient {exact} not in [{low},{high}]')


def variables(mx,my,hx=0.0,hy=0.0):
    return (Jet.variable(IA([[mx-hx]],[[mx+hx]]),0),
            Jet.variable(IA([[my-hy]],[[my+hy]]),1))


def power(jet,n):
    value=Jet.constant(1)
    for _ in range(n):value=value*jet
    return value


def main():
    Jet.configure(6)
    checked=0
    x,y=variables(1/8,1/16)
    inverse=(1-x-y).reciprocal()
    for (a,b),value in zip(Jet.indices,inverse.c):
        exact=F(math.comb(a+b,a))/F(13,16)**(a+b+1)
        enclosed(value,exact);checked+=1
    x,y=variables(1/8,1/16,1/64,1/128)
    inverse=(1-x-y).reciprocal()
    for (a,b),value in zip(Jet.indices,inverse.c):
        for base in (F(13,16)-F(3,128),F(13,16)+F(3,128)):
            enclosed(value,F(math.comb(a+b,a))/base**(a+b+1));checked+=1
    x,y=variables(1/4,5/16)
    root=(1+x+y).sqrt()
    for (a,b),value in zip(Jet.indices,root.c):
        n=a+b
        binomial=F(1)
        for k in range(n):binomial*= (F(1,2)-k)/(k+1)
        enclosed(value,binomial*math.comb(n,a)*F(5,4)*F(25,16)**(-n));checked+=1
    x,y=variables(0.,0.)
    arcsine=(x+2*y).asin()
    for (a,b),value in zip(Jet.indices,arcsine.c):
        n=a+b
        exact=F(0) if n%2==0 else F(math.comb(n-1,(n-1)//2),2**(n-1)*n)*math.comb(n,a)*2**b
        enclosed(value,exact);checked+=1

    mx,my,hx,hy=F(1,4),F(1,8),F(1,16),F(1,32)
    x,y=variables(float(mx),float(my))
    centers=power(1+x-2*y,8).lane_sum().c
    x,y=variables(float(mx),float(my),float(hx),float(hy))
    box=power(1+x-2*y,8).lane_sum().c
    integral=cell_enclosure(centers,box,0,float(hx),float(hy))
    exact=F(0)
    c=1+mx-2*my
    for a in range(9):
        for b in range(9-a):
            if a%2 or b%2:continue
            coefficient=F(math.factorial(8),math.factorial(a)*math.factorial(b)*math.factorial(8-a-b))*c**(8-a-b)*(-2)**b
            exact+=coefficient*4*hx**(a+1)*hy**(b+1)/((a+1)*(b+1))
    enclosed(integral,exact)
    print(json.dumps({'checks_passed':True,'exact_coefficient_endpoint_checks':checked,
                      'mixed_Taylor_remainder_polynomial_degree':8,
                      'exact_polynomial_integral':str(exact)},indent=2))


if __name__=='__main__':main()
