# Rigorous local certificates for the P4 Gaussian functional

## Certified results and their scope

This directory certifies the four-vertex path cost at two explicitly specified
points of the 4,655-type cone. It also certifies positive-radius neighborhoods
of those points. **It does not bound the entire cone.**

The nonuniform point is obtained by interpreting the six supported masses in
`../cone_search_P4.json` as exact decimal rationals, then normalizing their sum
exactly. The natural point uses the exact tree masses
`1/24, 1/24, 1/12, 1/6, 1/3, 1/3`, attached to the corresponding types.
The complete rational masses and types are saved in each certificate.

The computed per-vertex enclosures, rounded outward here, are

| Point | Lower endpoint | Upper endpoint |
|---|---:|---:|
| Rationalized nonuniform candidate | 0.136135422304858 | 0.136135429221129 |
| Exact natural tree | 0.136126085221635 | 0.136126385578789 |

In particular, the candidate exceeds the natural tree by **strictly more than
`9/1000000` per vertex**. Thus the natural tree is not a maximizer of this
relaxed P4 cone. No claim that the candidate itself is a maximizer is made.

Write `C_*` for either certified covariance matrix. Every symmetric matrix
with diagonal one and

\[
\max_{i<j}|C_{ij}-(C_*)_{ij}|\le2^{-20}
\]

is positive definite and satisfies

\[
\frac14\Phi_{P_4}(C)<\frac{6807}{50000}=0.13614.
\]

For the candidate, the verified intermediate upper bound is less than
`0.136139686695756`. Moreover, every mass vector in the 4,655-dimensional
simplex with

\[
\|m-m_*\|_1\le2^{-30}
\]

maps inside the certified covariance box and therefore satisfies the same
upper bound `0.13614`. This includes small off-support perturbations. The
radius and its normalization estimate are explicit; this remains a local
result.

The natural P4 point also has total expected cost less than `2 M_tree` by
strictly more than `433/50000 = 0.00866`, using the independently certified
tree constant in `../../mtree_refinement/certificate.json`.

All these comparisons are checked with exact rational arithmetic by
`compare_certificates.py`. The JSON files include exact rational versions of
the binary64 enclosure endpoints. Their outward decimal strings can be read
as rational decimal bounds; other JSON floating values should be interpreted
as binary64 diagnostics.

## Wider natural region used by the improved stability proof

[WIDE_NEIGHBORHOOD.md](WIDE_NEIGHBORHOOD.md) and
[wide_neighborhood.json](wide_neighborhood.json) certify a different,
weaker upper bound on a much wider region. Let `C_T` be the natural
boundary covariance and set `L=C_T-(J-I)/300`. Interval integration proves

```
0.553011755222534 <= Phi_4(L) <= 0.553012975913878
2*M_tree - Phi_4(L) > 1/10000.
```

Submodular Gaussian comparison extends the upper bound to every positive
semidefinite unit-diagonal `C` satisfying `C_ij >= (C_T)_ij-1/300`.
The upper coordinates are unrestricted apart from validity as a correlation
matrix. The full two-sided radius-`1/300` box is also certified positive
definite, but the theorem's comparison region is one-sided.

[MASS_NEIGHBORHOOD.md](MASS_NEIGHBORHOOD.md) and
[mass_neighborhood.json](mass_neighborhood.json) prove that every normalized
nonnegative mass vector within L1 distance `9/2500=0.0036` of the natural
mass vector maps into this comparison region. The interval derivative
calculation includes both normalizations and directions toward all 4,655
types. Its certified maximum downward correlation loss is less than
`0.003263615<1/300`.

The wider region supplies the local input to
[the improved global stability proof](../STABILITY_REFINEMENT.md).
It certifies `Phi_4<2*M_tree-1/10000`; it does not certify the stronger
per-vertex `0.13614` bound throughout the larger mass ball.

Replay from the repository root:

```sh
python -B research/gadget_insertion/certification/wide_neighborhood.py \
  --output research/gadget_insertion/certification/wide_neighborhood.json
python -B research/gadget_insertion/certification/mass_neighborhood.py \
  --output research/gadget_insertion/certification/mass_neighborhood.json
```

## Separate certificate for the H7 obstruction

[window_obstruction.py](window_obstruction.py) and
[window_obstruction.json](window_obstruction.json) certify a different
claim: an explicit positive rational point on the projected P8 tree face
has forward H7, reversed H7, and orientation-averaged expectation each
strictly greater than `17/125=0.136`. This exceeds the updated comparison
target `0.1341` and rules out the corresponding bare uniform window bound
on that relaxation.

The proof uses exact Boolean lower bounds by three-score events, both
covariance normalizations, interval Cholesky, and interval evaluation of
the Gaussian pair-sign arcsine formula. It requires no high-dimensional
Gaussian integration. It supplies no finite graph realization and does
not obstruct a bound averaged over consistent neighboring windows.
The complete argument and rational witness are in
[WINDOW_OBSTRUCTION.md](../WINDOW_OBSTRUCTION.md).

```sh
python -B research/gadget_insertion/certification/window_obstruction.py \
  --output research/gadget_insertion/certification/window_obstruction.json
```

## Reproduction of the original P4 certificates

Requirements are Python 3, NumPy, and mpmath. No SciPy integration routine or
LP solver is used. From the repository root:

```sh
python -B research/gadget_insertion/certification/verify_taylor.py

python -B research/gadget_insertion/certification/p4_interval.py \
  --degree 6 --tolerance 0.0000002 \
  --output research/gadget_insertion/certification/p4_point.json

python -B research/gadget_insertion/certification/p4_interval.py \
  --tree --degree 4 --tolerance 0.000004 \
  --output research/gadget_insertion/certification/p4_tree.json

python -B research/gadget_insertion/certification/compare_certificates.py \
  --output research/gadget_insertion/certification/comparison.json
```

The first successful complete runs took about 63 seconds and 52 seconds,
respectively. The candidate run closed 443 rectangles; the natural-tree run
closed 707. The exact sum of accepted rectangle areas was one in both runs.
Timing is machine-dependent. The tolerance controls subdivision, not the
validity of an enclosure: the actual returned interval is the certificate.

## 1. Exact Boolean data

For each of the 64 assignments to the six boundary occurrences, the program
enumerates all eight assignments to the three free path edges and counts
split vertices exactly. Its minimum is the integer cost `F(s)`.

For each pair of boundary coordinates `i,j`, define the mixed difference

\[
\Delta_{ij}F(s)=F(s+e_i+e_j)-F(s+e_i)-F(s+e_j)+F(s).
\]

The program evaluates all its 16 values on the remaining four coordinates.
It verifies that these differences are nonpositive and have absolute value
at most two. Its Walsh numerators are exact integers:

\[
n_A=\sum_{s\in\{0,1\}^4}\Delta_{ij}F(s)(-1)^{\sum_{a\in A}s_a}.
\]

Complement symmetry makes the odd-cardinality numerators zero. Therefore the
conditional expectation of a mixed difference requires only constants,
two-sign moments, and one four-sign moment. No numerical approximation is
used to obtain these coefficients or the independent value
`Phi(I_6)=sum_s F(s)/64`.

## 2. Reduction to a finite two-dimensional integral

For a positive-definite correlation matrix `C`, use the path

\[
C(t)=I+t(C-I),\qquad0\le t\le1.
\]

It remains positive definite. Gaussian differentiation and integration by
parts give the Plackett identity

\[
\frac{\partial\Phi}{\partial r_{ij}}
=\frac{\mathbb E[\Delta_{ij}F(Y)]}{2\pi\sqrt{1-r_{ij}^2}},
\]

where `Y` is the remaining Gaussian vector conditioned on the two selected
coordinates being zero. Signs are unchanged by normalizing the conditional
variances. Consequently

\[
\Phi(C)=\Phi(I_6)+\int_0^1\sum_{i<j}C_{ij}
 \frac{\partial\Phi}{\partial r_{ij}}(C(t))\,dt.
\]

The conditional four-dimensional correlations are obtained by Schur
complements. For conditioned pair `i,j` and remaining coordinates `k,l`, put
`r=C_ij`, with all correlations in this paragraph evaluated on the path.
The conditional correlation is

\[
\eta_{kl}=
\frac{C_{kl}(1-r^2)-C_{ki}C_{li}-C_{kj}C_{lj}
      +r(C_{ki}C_{lj}+C_{kj}C_{li})}
 {\sqrt{D_kD_l}},
\]

\[
D_k=1-r^2-C_{ki}^2-C_{kj}^2+2rC_{ki}C_{kj},
\]

and similarly for `D_l`. The numerators and these principal minors are
polynomials of degree at most three in `t`. The program expands them at the
cell midpoint and shifts the polynomial exactly to enclose the derivatives
throughout the cell. This retains cancellations without changing the formula.

For a four-dimensional correlation matrix `K`, let

\[
T_4(K)=\mathbb E\prod_{i=0}^3\operatorname{sign}(X_i).
\]

Scale only the four cross correlations between pairs `(0,1)` and `(2,3)` by
`u`. This is a convex path from the positive-definite block-diagonal matrix
to `K`, so it is positive definite throughout. Another application of the
same Gaussian differentiation identity gives

\[
T_4(K)=\frac4{\pi^2}\left[
 \arcsin K_{01}\arcsin K_{23}
 +\int_0^1\sum_{ij\in\{02,03,12,13\}}
   \frac{K_{ij}\arcsin\rho_{ij}(u)}{\sqrt{1-u^2K_{ij}^2}}\,du\right],
\]

where `rho_ij(u)` is the complementary conditional correlation.
The other sign moments are `(2/pi)*asin(K_ij)`.

Substituting the exact Walsh expansion into the outer derivative expresses
`Phi(C)` as its exact independent value plus a smooth integral over
`[0,1]^2`. Terms independent of `u` are integrated over a unit-length
interval. Both paths stay in the positive-definite domain, including at
their endpoints. Thus the integrand is smooth on a neighborhood of the
compact square. There is no truncated Gaussian tail or unproved quadrature
error in this reduction.

The code uses expanded inner-minor polynomials to preserve cancellations.
If the unscaled inner correlations are denoted by
`r=K_ij`, `z=K_kl`, `p=K_ki`, `q=K_lj`, `v=K_kj`, `w=K_li`, then

\[
\begin{aligned}
N(u)&=u(z-pw-vq+rpq)+u^3(rvw-zr^2),\\
D_k(u)&=1-p^2-u^2(r^2+v^2-2rpv),\\
D_l(u)&=1-q^2-u^2(r^2+w^2-2rqw).
\end{aligned}
\]

These are the same Schur-complement expressions, after collecting powers.

## 3. Interval Taylor coefficients and explicit remainders

A jet stores all normalized mixed derivatives

\[
f_{ab}=\frac{\partial_t^a\partial_u^bf}{a!b!},\qquad a+b\le N.
\]

Addition and multiplication use the usual convolution. Reciprocals and
square roots are computed by equating coefficients in `fg=1` and `g^2=f`.
Arcsine coefficients are obtained by differentiating
`asin(f)`, using `f'/sqrt(1-f^2)`. These recurrences hold at every base point;
performing them in interval arithmetic therefore encloses the coefficients
at every point of a rectangle.

For tighter bounds the program intersects intermediate coefficients with
their multivariate mean-value forms. If `c` is the rectangle center and
`delta` ranges over its centered rectangle, then

\[
f_{ab}(c+\delta)\in f_{ab}(c)
 +(a+1)[f_{a+1,b}]\delta_t
 +(b+1)[f_{a,b+1}]\delta_u.
\]

The brackets enclose the indicated coefficients throughout the rectangle.
Descending through derivative orders permits already tightened higher
derivatives to tighten lower ones. Intersecting two containing intervals
preserves containment. This step is essential for computational efficiency,
but not a heuristic error estimate.

For an even order `N` and a rectangle centered at `c` with halfwidths
`h_t,h_u`, the multivariate Taylor theorem gives

\[
f(c+\delta)=\sum_{a+b<N}f_{ab}(c)\delta_t^a\delta_u^b
 +\sum_{a+b=N}f_{ab}(c+\theta\delta)\delta_t^a\delta_u^b,
\]

for some `0<=theta<=1` depending on `delta`. Integrating the terms of degree
less than `N` is exact up to containing arithmetic: an odd exponent gives
zero, while an even-even monomial has integral

\[
J_{ab}=\frac{4h_t^{a+1}h_u^{b+1}}{(a+1)(b+1)}.
\]

For a remainder coefficient known to lie in `[L_ab,U_ab]` throughout the
rectangle:

* If both exponents are even, its integrated contribution lies in
  `[L_ab,U_ab]*J_ab` because the monomial is nonnegative.
* If both exponents are odd, its contribution lies in
  `[-(U_ab-L_ab)*J_ab/2, (U_ab-L_ab)*J_ab/2]`. Subtract the midpoint of the
  coefficient interval; its constant part integrates to zero, and the
  remaining absolute value is at most half the interval width.

Since `N` is even, these are the only possibilities in the remainder.
The program does **not** discard odd remainder terms merely because their
monomials integrate to zero with a constant coefficient.

The initial rectangles tile the unit square. A rejected rectangle is replaced
by its two exact dyadic halves, preserving a disjoint cover. Only rectangles
with finite containing bounds are accepted. The algorithm terminates
successfully only when no rectangle remains, and also checks that their
accepted areas sum to one as an exact rational number. Adding their interval
integrals therefore gives a containing enclosure of the whole integral.

## 4. Positive definiteness and the local correlation box

The input masses are exact positive rationals. Powers of `sqrt(2)`, their Gram
matrix, and both covariance normalizations are enclosed by interval
arithmetic. Interval Cholesky with strictly positive diagonal pivots proves
that every matrix in the resulting entry box is positive definite.

Let `rho=2^-20`. The program also performs interval Cholesky on
`C_* - 5*rho*I`. For a symmetric perturbation with zero diagonal and each
off-diagonal entry of absolute value at most `rho`, its operator norm is at
most `5*rho` by the maximum absolute row-sum bound. The entire specified
correlation box is therefore positive definite.

Exact Boolean enumeration gives `|Delta_ij F|<=2`. Throughout this box the
program checks `|r_ij|<24/25`. The derivative formula implies

\[
\left|\frac{\partial\Phi}{\partial r_{ij}}\right|
\le\frac1{\pi\sqrt{1-r_{ij}^2}}
\le\frac{25}{7\pi}<\frac{25}{21}.
\]

There are 15 upper-triangle entries. Dividing the integrated derivative bound
by four vertices gives

\[
\frac14\Phi(C)\le\frac14\Phi(C_*)+\frac{125}{28}\rho.
\]

The rational comparison of the computed upper endpoint with `6807/50000` is
checked explicitly. This proves the advertised upper bound throughout the
box, not merely at the center.

## 5. A local neighborhood in the full mass simplex

Every generator of `T_4` has coordinates `g_t(i)=2^{-t_i/2}` with
`|t_i|<=4`, by the required-edge distance inequalities. Hence
`0<g_t(i)g_t(j)<=16`. For normalized nonnegative masses,

\[
\|m-m_*\|_1\le\delta\quad\Longrightarrow\quad
\max_{ij}|G(m)_{ij}-G(m_*)_{ij}|\le16\delta.
\]

The center's Gram diagonals are exact rationals
`sum_t m_t*2^(-t_i)`. The program computes them exactly and checks that every
one is at least `4/5`.

The normalization estimate used twice is as follows. Suppose `A` and `B` are
positive-semidefinite matrices, `min A_ii>=d`, and
`max|B_ij-A_ij|<=eta`. Put `r=eta/d<=1/2`. Then all diagonals of `B` are
positive, and

\[
\max_{ij}|\operatorname{corr}(B)_{ij}-\operatorname{corr}(A)_{ij}|
\le\frac{2r}{1-r}\le4r.
\]

Indeed, write the normalized perturbed entry as
`(corr(A)_ij+e_ij)/sqrt((1+e_i)(1+e_j))`, with every error bounded by `r`.
Use `|corr(A)_ij|<=1` and the bounds `[1/(1+r),1/(1-r)]` on the reciprocal
square root.

The first normalization therefore changes vertex correlations by at most
`80*delta`. Each edge-key Gram entry is a sum of four such entries and changes
by at most `320*delta`. Its center diagonal is at least two because all
generator coordinates are positive. The second normalization changes key
correlations by at most `640*delta`, provided the indicated relative errors
are at most one half.

Taking `delta=2^-30` satisfies these conditions and
`640*delta<2^-20`. Thus the entire declared L1 neighborhood in the full
simplex maps into the certified correlation box. No support restriction is
imposed on the nearby mass vector.

## 6. Trust base and verification

The interval arithmetic is the pinned upstream `mtree/enclosure/ivarith.py`,
Git blob `b99d01f757f2b3d75151ba9403a547c6b8c60477`. The program verifies that
hash. Each binary64 arithmetic operation is followed by one outward ULP;
square roots use NumPy's square root with outward rounding. Arcsine endpoints
are checked against 80-bit mpmath interval sine, rather than trusting an
ordinary inverse-trigonometric evaluation. The enclosure of pi is also
checked against mpmath's interval constant.

The arithmetic assumptions are the same IEEE round-to-nearest and NumPy
binary64 assumptions documented in the upstream module. All structural,
positivity, scope, and target comparisons raise explicit exceptions rather
than relying on assertions that `python -O` could remove.

`verify_taylor.py` checks 140 Taylor coefficients or coefficient endpoints
against exact rational formulas for a reciprocal, a square root at a rational
square, and an arcsine at zero. It also encloses an independently integrated
degree-eight polynomial using the same mixed Taylor remainder routine.
Separate numerical spot checks agreed with the independent
`../path_gradient.py` derivative evaluator to approximately `5e-16`; these
spot checks are diagnostics and are not part of the proof.

These computations do not prove global optimality, do not cover the remainder
of the P4 cone, and do not certify P5, P6, or P7. Their proved contributions
are the explicit point enclosures, strict point comparisons, and the stated
local neighborhoods.
