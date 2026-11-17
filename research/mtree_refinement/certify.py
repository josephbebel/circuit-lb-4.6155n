#!/usr/bin/env python3
"""Certify the tree constant with integer interval arithmetic only.

No third-party libraries, numerical quadrature, or floating-point arithmetic are
used in the certificate.  See README.md in this directory for the identities,
the positive-series construction, and the explicit remainder estimate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from math import isqrt
from pathlib import Path


SCALE = 10**70


def ceildiv(n: int, d: int) -> int:
    return -((-n) // d)


@dataclass(frozen=True)
class I:
    """The closed rational interval [lo/SCALE, hi/SCALE]."""

    lo: int
    hi: int

    def __post_init__(self):
        if self.lo > self.hi:
            raise ArithmeticError("reversed interval")

    @staticmethod
    def of(value: I | int) -> I:
        if isinstance(value, I):
            return value
        if not isinstance(value, int):
            raise TypeError("only exact integers or intervals are allowed")
        return I(value * SCALE, value * SCALE)

    def __add__(self, other):
        other = I.of(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-I.of(other))

    def __rsub__(self, other):
        return I.of(other) + (-self)

    def __mul__(self, other):
        other = I.of(other)
        products = (
            self.lo * other.lo, self.lo * other.hi,
            self.hi * other.lo, self.hi * other.hi,
        )
        return I(min(products) // SCALE, ceildiv(max(products), SCALE))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ArithmeticError("division by an interval containing zero")
        return I(SCALE**2 // self.hi, ceildiv(SCALE**2, self.lo))

    def __truediv__(self, other):
        return self * I.of(other).reciprocal()

    def __rtruediv__(self, other):
        return I.of(other) * self.reciprocal()

    def sqrt(self):
        if self.lo < 0:
            raise ArithmeticError("square root of a negative interval")
        low = isqrt(self.lo * SCALE)
        high = isqrt(self.hi * SCALE)
        if high * high != self.hi * SCALE:
            high += 1
        return I(low, high)

    def __pow__(self, exponent: int):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("only nonnegative integer powers are allowed")
        result = I.of(1)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent //= 2
        return result


def require(condition: bool, message: str):
    if not condition:
        raise ArithmeticError(message)


def rational(n: int, d: int = 1) -> I:
    require(d > 0, "rational denominator must be positive")
    return I((n * SCALE) // d, ceildiv(n * SCALE, d))


def arctan_reciprocal(q: int, terms: int) -> I:
    """Alternating arctan(1/q) series with its next-term remainder."""
    require(q > 1 and terms > 0, "bad arctangent parameters")
    total = I.of(0)
    power = q
    for n in range(terms):
        term = rational(1, (2*n + 1) * power)
        total = total + term if n % 2 == 0 else total - term
        power *= q*q
    error = rational(1, (2*terms + 1) * power)
    # The true remainder has the sign of the first omitted term.
    return total + (I(0, error.hi) if terms % 2 == 0 else I(-error.hi, 0))


def pi_interval(terms: int) -> I:
    # Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    return 16*arctan_reciprocal(5, terms) - 4*arctan_reciprocal(239, terms)


def asin_positive(x: I, terms: int) -> I:
    """Positive arcsine series, with a geometric majorant for its tail."""
    require(x.lo >= 0 and x.hi < SCALE, "arcsine argument not in [0,1)")
    require(terms > 0, "arcsine requires at least one term")
    xx = x*x
    central = I.of(1)
    power = x
    total = I.of(0)
    for n in range(terms):
        total += central*power/(2*n + 1)
        central = central*(2*n + 1)/(2*n + 2)
        power *= xx
    # c_n = binom(2n,n)/4^n <= 1, hence the omitted positive terms
    # are at most x^(2*terms+1)/(1-x^2).
    tail = power/(1-xx)
    return I(total.lo, total.hi + tail.hi)


def asin_large(x: I, pi: I, terms: int) -> I:
    require(x.lo >= 0 and x.hi < SCALE, "half-angle input not in [0,1)")
    return pi/2 - 2*asin_positive(((1-x)/2).sqrt(), terms)


def positive_integral(p: I, s: I, terms: int):
    """Enclose integral_0^s asin(aq/(1-bq^2))/sqrt(1-q^2) dq."""
    require(terms > 0, "integral requires at least one term")
    a = (1-p)/(1+p)
    b = 2/(1+p)
    bb = b*b
    require(a.lo > 0 and b.lo > SCALE, "positive-series parameter failure")

    # c_n are the central binomial coefficients divided by 4^n.
    central = [I.of(1)]
    weighted = [I.of(1)]
    bpower = I.of(1)
    for n in range(1, terms):
        central.append(central[-1]*(2*n-1)/(2*n))
        bpower *= bb
        weighted.append(central[-1]*bpower)

    # e_n are coefficients of (1-x)^(-1/2)(1-b^2 x)^(-1/2).
    # d_n are coefficients of y'(q), x=q^2,
    # y(q)=asin(aq/(1-bq^2)).
    derivative = []
    geometric = I.of(0)
    previous_e = I.of(0)
    for n in range(terms):
        e = sum((central[j]*weighted[n-j] for j in range(n+1)), I.of(0))
        if n:
            geometric = b*(previous_e + geometric)
        derivative.append(a*(e + 2*geometric)/(2*n+1))
        previous_e = e

    # derivative[n] now represents the coefficient of q^(2n+1) in y.
    ss = s*s
    power = ss
    partial = I.of(0)
    for n in range(terms):
        coefficient = sum(
            (derivative[j]*central[n-j] for j in range(n+1)), I.of(0)
        )
        require(coefficient.lo >= 0, "coefficient failed positivity check")
        partial += coefficient*power/(2*n+2)
        power *= ss

    # h(q)/q = sum c_n q^(2n) with c_n >= 0.  Bound it at R=39/40.
    radius = rational(39, 40)
    require(0 < s.lo and s.hi < radius.lo, "tail radius does not exceed s")
    require((b*radius).hi < SCALE, "tail radius crosses a singularity")
    rho = a*radius/(1-b*radius*radius)
    require(0 <= rho.lo and rho.hi < SCALE, "tail arcsine argument invalid")
    # asin(rho) <= rho/sqrt(1-rho^2), by integrating an increasing derivative.
    h_over_radius = rho/(radius*(1-rho*rho).sqrt()*(1-radius*radius).sqrt())
    require(h_over_radius.hi < 10*SCALE, "tail majorant 10 was not certified")
    tail = 10*ss*((s/radius)**(2*terms))/(2*terms+2)
    result = I(partial.lo, partial.hi+tail.hi)
    return result, tail, h_over_radius


def fixed_decimal(n: int, digits: int) -> str:
    sign = "-" if n < 0 else ""
    n = abs(n)
    unit = 10**digits
    return f"{sign}{n//unit}.{n%unit:0{digits}d}"


def outward_decimal(x: I, scale_digits: int, output_digits: int) -> list[str]:
    require(output_digits <= scale_digits, "more output digits than working precision")
    divisor = 10**(scale_digits-output_digits)
    return [
        fixed_decimal(x.lo//divisor, output_digits),
        fixed_decimal(ceildiv(x.hi, divisor), output_digits),
    ]


def interval_record(x: I, digits: int, output_digits: int) -> dict:
    return {
        "decimal_enclosure": outward_decimal(x, digits, output_digits),
        "lower_numerator": str(x.lo),
        "upper_numerator": str(x.hi),
        "denominator": str(SCALE),
    }


def certify(digits: int = 70, terms: int = 320, output_digits: int = 30) -> dict:
    global SCALE
    require(digits >= 40 and terms >= 1 and output_digits >= 1, "invalid precision")
    SCALE = 10**digits
    transcendental_terms = digits + 20
    pi = pi_interval(transcendental_terms)
    require(3*SCALE < pi.lo < pi.hi < 4*SCALE, "pi enclosure failed")
    root2 = I.of(2).sqrt()
    p = (1+2*root2)/4
    s = (2+root2)/4
    asin_p = asin_large(p, pi, transcendental_terms)
    asin_s = asin_large(s, pi, transcendental_terms)
    correction, tail, majorant = positive_integral(p, s, terms)
    independent_pairs_value = 1-2*asin_p/pi+2*(rational(1,4)+asin_p/(2*pi))**2
    mtree = independent_pairs_value-asin_s/pi+2*correction/(pi*pi)
    require(mtree.lo > 0, "nonpositive M_tree enclosure")
    kappa = mtree/2
    circuit = 1+1/mtree
    old_upper = rational(27658360966404055, 10**17)
    require(mtree.hi < old_upper.lo, "did not improve the old certified upper bound")
    return {
        "certificate": "M_tree positive power series with exact integer interval arithmetic",
        "schema_version": 1,
        "upstream_commit": "a56500014ba159397903dea81d41ac836f6880d5",
        "working_decimal_places": digits,
        "integral_terms": terms,
        "transcendental_terms": transcendental_terms,
        "M_tree": interval_record(mtree, digits, output_digits),
        "cubic_cutwidth_coefficient_M_tree_over_2": interval_record(kappa, digits, output_digits),
        "circuit_coefficient_1_plus_1_over_M_tree": interval_record(circuit, digits, output_digits),
        "integral_tail_bound": interval_record(I(0, tail.hi), digits, output_digits),
        "tail_majorant_bound": interval_record(majorant, digits, output_digits),
        "claims": {
            "M_tree_upper_improves_original": True,
            "changes_matched_insertion_method": False,
            "new_leaderboard_record": False,
            "full_upstream_certificates_replayed": False,
            "circuit_and_cutwidth_consequences_depend_on_upstream_theorems": True,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=70)
    parser.add_argument("--terms", type=int, default=320)
    parser.add_argument("--output-digits", type=int, default=30)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = certify(args.digits, args.terms, args.output_digits)
    text = json.dumps(record, indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
