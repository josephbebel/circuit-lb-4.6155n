# A wide certified neighborhood of the natural P4 point

## Result

Let `C_T` be the exact six-by-six boundary-score correlation matrix of the
natural P4 distance-type masses, and let `Phi_4` denote the **total** optimized
four-vertex insertion cost. Put

\[
r=\frac1{300},\qquad L=C_T-r(J-I).
\]

The interval calculation in `wide_neighborhood.py` proves

\[
0.553011755222534
\le\Phi_4(L)\le
0.553012975913878,
\]

where these decimal endpoints are rounded outward. The JSON result stores
the exact rational versions of the binary64 endpoints.

Consequently, for **every** positive-semidefinite correlation matrix `C`
with unit diagonal satisfying

\[
C_{ij}\ge (C_T)_{ij}-\frac1{300}\quad(i<j),                 \tag{1}
\]

we have

\[
\Phi_4(C)\le\Phi_4(L)<2M_{\rm tree}-\frac1{10000}.         \tag{2}
\]

This is a one-sided region; upper entry deviations are unrestricted except
for the correlation-matrix conditions. In particular, it contains the full
two-sided entrywise covariance box of radius `1/300` around `C_T`. Every
symmetric unit-diagonal matrix in that two-sided box is positive definite,
as certified below.

It also contains the image of the **full normalized mass-simplex ball**

\[
\|m-m_{\rm tree}\|_1\le\frac1{10000}.                    \tag{3}
\]

This is over 100,000 times the former certified mass radius `2^-30`.
Its cost bound is deliberately weaker than the earlier `0.13614` per-vertex
target: the purpose is to certify a much larger region with a positive
deficit relative to two old pair gadgets. It does not bound the entire
4,655-type simplex.

## 1. Why one comparison point suffices

The Boolean P4 cost is submodular. For each pair of boundary coordinates,
the program enumerates the mixed difference

\[
\Delta_{ij}F(s)=F(s+e_i+e_j)-F(s+e_i)-F(s+e_j)+F(s)
\]

on all 16 assignments to the remaining four coordinates and checks that
every value is nonpositive. Gaussian differentiation therefore gives

\[
\frac{\partial\Phi_4}{\partial C_{ij}}
=\frac{\mathbb E[\Delta_{ij}F(Y)]}
       {2\pi\sqrt{1-C_{ij}^2}}\le0
\]

at every positive-definite correlation matrix. Here the expectation is
conditional on the selected two Gaussian coordinates being zero.

The program proves `L` positive definite by interval Cholesky. If `C`
satisfies (1), the segment `(1-t)L+tC` consists of positive-definite
correlation matrices for `0<=t<1`, and every off-diagonal entry is
nondecreasing along it. Thus the cost is nonincreasing. Continuity of
Gaussian sign expectations at unit-diagonal positive-semidefinite limits
extends the inequality to a possibly singular endpoint `C`.

The numerical enclosure is only needed at `L`; it is not necessary to
integrate an entire interval covariance box. This monotonic comparison
avoids the losses in the earlier uniform derivative estimate.

## 2. The rigorous point calculation and positive definiteness

The center is constructed from the exact natural masses

\[
(1/24,1/24,1/12,1/6,1/3,1/3)
\]

on the six specified tree types. Their total is one and all ten Gram
diagonals are exactly one; the program verifies those rational identities.
Both normalizations and subtraction of the exact rational `1/300` are
performed with containing interval arithmetic.

The evaluator is the previously audited `p4_interval.py`: an exact Boolean
Walsh expansion followed by a finite two-dimensional Plackett integral and
explicit multivariate Taylor remainder bounds. Its proof, including the
odd-odd top-degree remainder terms and mean-value tightening, is in
`README.md`, Sections 1-3. This calculation uses degree four and an adaptive
dyadic cover. The successful run closed 645 rectangles, with exact total
accepted area one. The stopping tolerance controls subdivision; validity
comes from the returned containing intervals.

The code additionally performs interval Cholesky on

\[
C_T-5rI.
\]

An arbitrary symmetric zero-diagonal perturbation whose entries have
absolute value at most `r` has operator norm at most `5r`. Positive
definiteness of the shifted center therefore proves positive definiteness
throughout the complete two-sided covariance box.

The lower enclosure of `M_tree` from
`../../mtree_refinement/certificate.json` gives, using exact fractions,

\[
2M_{\rm tree}-\Phi_4(L)>\frac1{10000}.
\]

The exact difference is recorded in `wide_neighborhood.json`.

## 3. A larger full-simplex mass neighborhood

For a normalized P4 mass vector, write

\[
G(m)_{ij}=\sum_t m_t g_t(i)g_t(j),\qquad
g_t(i)=2^{-t_i/2}.
\]

Every generator coordinate satisfies `|t_i|<=4`, so every product lies in
`[1/16,16]`. The program checks this integer coordinate bound over all
4,655 types.

Let `d=||m-m_tree||_1`. Because both mass vectors sum to one, their signed
difference has equal positive and negative total masses `d/2`. A linear
functional whose coefficients lie in `[a,b]` changes by at most
`(b-a)d/2`. Therefore

\[
\max_{ij}|G(m)_{ij}-G(m_{\rm tree})_{ij}|
\le\frac{255}{32}d<8d.                                  \tag{4}
\]

The zero-sum condition improves the earlier bound `16d` without assuming
that the perturbation remains on the six-type face.

For completeness, if two positive-semidefinite matrices have center
diagonals at least `D` and entrywise difference at most `e<D`, their
normalized correlations differ by at most

\[
\frac{2e}{D-e}.                                         \tag{5}
\]

Indeed, divide the perturbed entry and its two diagonals by the center
standard deviations. The numerator changes by at most `e/D`, and the
reciprocal normalization factor lies between `1/(1+e/D)` and
`1/(1-e/D)`. The center normalized entry has absolute value at most one,
which gives (5).

At the natural P4 point all vertex Gram diagonals equal one. Put
`epsilon=8d`. The first normalization changes every vertex correlation by
at most

\[
v=\frac{2\epsilon}{1-\epsilon}.
\]

Each boundary-key Gram entry is the sum of four vertex correlations, so
its change is at most `E=4v`. Its natural diagonal is

\[
D=2+2\rho=2+\frac{4\sqrt2}{3}>\frac{97}{25}=3.88.
\]

The last inequality follows from `sqrt(2)>141/100`, checked by squaring
positive rationals. A second application of (5) bounds the key-correlation
change by

\[
\frac{2E}{97/25-E}.
\]

For `d<=1/10000`, exact arithmetic gives

\[
\epsilon\le\frac1{1250},\qquad
v\le\frac2{1249},\qquad
E\le\frac8{1249},\qquad
\max_{ij}|C_4(m)_{ij}-(C_T)_{ij}|
\le\frac{400}{120953}<\frac1{300}.                       \tag{6}
\]

Thus (3) implies (1) and the cost deficit (2), including arbitrary
off-support masses allowed by the full simplex.

## 4. Replay and scope

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -B \
  research/gadget_insertion/certification/wide_neighborhood.py \
  --output research/gadget_insertion/certification/wide_neighborhood.json
```

Requirements and interval-arithmetic assumptions are the same as for the
existing P4 certificates: Python, NumPy, mpmath, IEEE binary64 arithmetic
with outward rounding, and verified inverse-trigonometric endpoints. The
script checks the pinned upstream interval-arithmetic Git blob and verifies
its pi interval against mpmath's interval constant. Certificate inputs are
identified by hashes in the result.

This result is a local theorem with an explicit, substantially enlarged
region. A universal conclusion requires a separate argument putting the
remaining cone points into this region or bounding their cost by another
method. No global optimizer or leaderboard claim follows from this file
alone.
