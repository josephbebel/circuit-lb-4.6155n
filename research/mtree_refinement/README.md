# A small exact-arithmetic certificate for the tree constant

This directory improves the numerical enclosure of the existing tree constant.
It does not change the matched-insertion ordering or its extremal value.
The verifier uses Python's standard library, integer arithmetic, and explicitly
bounded convergent series. It does not use floating-point quadrature, an LP
solver, or a transcendental-function library.

## Result and reproduction

Run from the repository root:

```sh
python -B research/mtree_refinement/certify.py \
  --output research/mtree_refinement/certificate.json
```

The default calculation takes approximately a fraction of a second in the
development environment. Its outward-rounded rational decimal enclosures are

\[
\begin{aligned}
0.276583606882373471022763934877
&\le M_{\rm tree}\le
0.276583606882373471022763934878,\\
0.138291803441186735511381967438
&\le M_{\rm tree}/2\le
0.138291803441186735511381967439,\\
4.615543275582792605180647436566
&\le 1+1/M_{\rm tree}\le
4.615543275582792605180647436567.
\end{aligned}
\]

The first enclosure replaces the numerical precision of the upstream interval
`[0.27658360687864314, 0.27658360966404055]`. The JSON also contains the exact
integer numerators and common denominator of the narrower internal interval.
Its width is less than `2.435e-40`; the displayed endpoints are rounded outward
to 30 places. The integral-tail upper bound is less than `1.202e-39`.

The relevant upstream version is
[`a56500014ba159397903dea81d41ac836f6880d5`](https://github.com/rrwilliams0/circuit-lb-4.6155n/tree/a56500014ba159397903dea81d41ac836f6880d5).
The definitions, Walsh identity, and original numerical enclosure are in
[`proof.md`, Sections 4.1–4.3](https://github.com/rrwilliams0/circuit-lb-4.6155n/blob/a56500014ba159397903dea81d41ac836f6880d5/proof.md).

### What this establishes

This is an independent certificate of the numerical value of the explicit
Gaussian functional at the tree matrix. Combined with the upstream theorem
`M_123 = M_tree` and its circuit-to-cutwidth reduction, it strengthens the
certified circuit coefficient from `4.6155432392` to at least
`4.615543275582792605180647436566`.

The circuit and graph-theoretic conclusions still depend on the upstream
theorems and their finite certificates. This calculation does not replay those
large certificates, certify a new ordering rule, or establish a new leaderboard
record. Upstream Section 1.2 proves that the present matched-insertion cone
cannot give a coefficient greater than the exact value `1 + 1/M_tree`.

No upstream source or old certificate needs to be replaced: the original
conservative numerical constants remain valid throughout the original proof.

## 1. An exact one-dimensional formula

Put

\[
p=\frac{1+2\sqrt2}{4},\qquad
s=\frac{2+\sqrt2}{4}.
\]

For `0 <= q <= s`, consider the correlation matrix

\[
C(q)=\begin{pmatrix}
1&p&q&q\\p&1&q&q\\q&q&1&p\\q&q&p&1
\end{pmatrix}.
\]

Its eigenvalues are `1-p`, `1-p`, `1+p-2q`, and `1+p+2q`.
All are positive on this interval; in particular `1+p-2s = 1/4`.
Let `F(q)` be the Gaussian matched-edge functional of upstream Section 4.2,
and let `Q(q)` be the positive-orthant probability for `N(0,C(q))`.
The Walsh identity gives

\[
F(q)=1-\frac{2}{\pi}\bigl(\arcsin p+\arcsin q\bigr)+2Q(q).
\]

At `q=0`, the two pairs are independent, so

\[
F(0)=1-\frac{2\arcsin p}{\pi}
 +2\left(\frac14+\frac{\arcsin p}{2\pi}\right)^2.
\]

For a cross pair, conditioning both of its Gaussian coordinates to zero gives
the remaining pair the correlation

\[
\rho(q)=\frac{q(1-p)}{1+p-2q^2}.
\]

This follows directly from the Schur complement of the covariance matrix of
the conditioned coordinates. Differentiating the Gaussian density with
respect to a correlation and integrating over the orthant gives the Plackett
identity: the derivative is the density of the conditioned pair at zero times
the conditional orthant probability. There are four equal cross derivatives.
The zero density is `1/(2*pi*sqrt(1-q^2))`; the conditional orthant probability
is `1/4 + asin(rho(q))/(2*pi)`. Therefore

\[
Q'(q)=\frac{1}{2\pi\sqrt{1-q^2}}
       +\frac{\arcsin\rho(q)}{\pi^2\sqrt{1-q^2}},
\]

and hence

\[
F'(q)=-\frac{1}{\pi\sqrt{1-q^2}}
       +\frac{2\arcsin\rho(q)}{\pi^2\sqrt{1-q^2}}.
\]

All differentiations occur inside the positive-definite domain; on this
compact path, standard Gaussian domination justifies integration and
differentiation. Integrating the last identity yields the formula certified by
the program:

\[
\boxed{
M_{\rm tree}=F(0)-\frac{\arcsin s}{\pi}
 +\frac{2}{\pi^2}\int_0^s
   \frac{\arcsin\!\left(aq/(1-bq^2)\right)}{\sqrt{1-q^2}}\,dq
}
\]

where

\[
a=\frac{1-p}{1+p},\qquad b=\frac{2}{1+p}=1+a.
\]

## 2. Positive coefficients

Let `y(q) = asin(a*q/(1-b*q^2))`. The algebraic identity

\[
(1-bq^2)^2-a^2q^2=(1-q^2)(1-b^2q^2)
\]

follows from `a=b-1`. For `0 <= q < 1/b`, it gives

\[
y'(q)=\frac{a(1+bq^2)}
 {(1-bq^2)\sqrt{1-q^2}\sqrt{1-b^2q^2}}.
\]

Every factor on the right has a power series in `q^2` with nonnegative
coefficients. Define

\[
c_n=\binom{2n}{n}/4^n,\qquad
c_0=1,\qquad c_{n+1}=c_n\frac{2n+1}{2n+2}.
\]

Writing `x=q^2`, the program first forms

\[
e_n=\sum_{j=0}^n c_jc_{n-j}b^{2(n-j)}.
\]

These are the coefficients of
`(1-x)^(-1/2)*(1-b^2*x)^(-1/2)`.
The coefficients of `(1+b*x)/(1-b*x)` are `1` at degree zero and `2*b^k`
at degree `k>=1`. Thus the coefficient of `q^(2n)` in `y'(q)` is

\[
d_n=a\left(e_n+2\sum_{k=1}^n b^k e_{n-k}\right).
\]

Integrating `y'` term by term and multiplying by `(1-q^2)^(-1/2)` shows

\[
h(q):=\frac{y(q)}{\sqrt{1-q^2}}
     =q\sum_{n=0}^\infty A_nq^{2n},\qquad
A_n=\sum_{j=0}^n\frac{d_jc_{n-j}}{2j+1}\ge0.
\]

All the series converge for `|q| < 1/b`, and `s < 1/b`.
Consequently the integral is

\[
I=\sum_{n=0}^\infty\frac{A_n s^{2n+2}}{2n+2}.
\]

The code computes the first `N=320` terms. The construction uses positive
convolutions and an accumulating geometric recurrence, avoiding cancellation
in the coefficients.

## 3. An explicit bound on the entire omitted tail

Set `R=39/40`. The verifier checks, using rational intervals, that

\[
0<s<R<1/b,\qquad 0\le\rho(R)<1.
\]

For `0 <= u < 1`, the increasing derivative of arcsine gives

\[
\arcsin u\le \frac{u}{\sqrt{1-u^2}}.
\]

It follows that

\[
\sum_{n=0}^\infty A_nR^{2n}
=\frac{h(R)}R
\le\frac{\rho(R)}{R\sqrt{1-\rho(R)^2}\sqrt{1-R^2}}<10.
\]

The verifier explicitly checks the last strict inequality. The expression
actually lies in `[5.213356417822692195899169859858,
5.213356417822692195899169859859]`.

Since all coefficients are nonnegative, the remainder after `N` terms obeys

\[
\begin{aligned}
0\le I-\sum_{n=0}^{N-1}\frac{A_n s^{2n+2}}{2n+2}
&\le\frac{s^2}{2N+2}
       \left(\frac{s}{R}\right)^{2N}
       \sum_{n=N}^\infty A_n R^{2n}\\
&\le\boxed{\frac{10s^2}{2N+2}
       \left(\frac{s}{R}\right)^{2N}}.
\end{aligned}
\]

This bounds the whole infinite tail. There is no empirical stopping criterion
or inference from the agreement of two numerical evaluations.

## 4. The arithmetic trust base

An interval is represented by two integers `L,U` and the common exact scale
`S=10^70`, meaning the rational interval `[L/S,U/S]`.

* Addition and negation are exact.
* Multiplication takes the minimum and maximum of the four integer endpoint
  products and rounds division by `S` down and up.
* Reciprocation rounds `S^2/U` down and `S^2/L` up after rejecting intervals
  containing zero.
* Square roots use `math.isqrt`: the lower endpoint is
  `floor(sqrt(L*S))`, and the upper endpoint is
  `ceil(sqrt(U*S))`. Thus both square-root inequalities are integer-checkable.
* Integer divisions and all reported decimal endpoints are rounded outward.
  Python's arbitrary-precision integers do not overflow.

Pi is enclosed by Machin's identity

\[
\pi=16\arctan(1/5)-4\arctan(1/239),
\]

using alternating arctangent series and the first omitted term as the error
bound. The arcsine values at `p` and `s` use

\[
\arcsin z=\frac\pi2-2\arcsin\sqrt{\frac{1-z}{2}},
\]

followed by the positive arcsine series

\[
\arcsin x=\sum_{n=0}^\infty\frac{c_n x^{2n+1}}{2n+1}.
\]

After `K` terms, its remainder is at most
`x^(2K+1)/(1-x^2)`, since `0 < c_n/(2n+1) <= 1`.
The default is `K=90`. Each power, coefficient, partial sum, and remainder is
itself computed as a containing rational interval.

Thus the numerical trust base is the readable verifier, Python's integer
operations and integer square root, and the analytic identities proved above.
The input constants are algebraic expressions, not quoted decimal values.

## 5. Validation and scope

The following checks were performed during development:

* Exact `fractions.Fraction` endpoint comparisons for interval addition,
  multiplication, and division, including positive, negative, mixed-sign,
  and zero intervals; exact square inequalities for square-root enclosures.
* A second rigorous run with 85 working decimal places and 400 integral
  terms. Its interval lies strictly inside the default interval and encloses
  `M_tree` between
  `0.2765836068823734710227639348778729937589` and
  `0.2765836068823734710227639348778729937590`.
* An independent 65-decimal numerical quadrature along the different path
  that scales all six off-diagonal correlations. It gives
  `0.276583606882373471022763934877872993758963969073626876185586`,
  which lies inside the integer certificate. This last check is a numerical
  sanity check and is not used as a proof step.

Additional digits can be requested, for example:

```sh
python -B research/mtree_refinement/certify.py \
  --digits 85 --terms 400 --output-digits 40
```

The program always returns containing intervals; increasing the number of
printed digits alone does not narrow them. All decisive validity checks raise
exceptions explicitly, so they remain active under `python -O`.
