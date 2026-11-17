# Sharper stability constants for the global P4 improvement

## Result and dependencies

This note strengthens the uniform P4 deficit from `10^-54` to

\[
M_4\le 2M_{\rm tree}-\frac{17}{5\cdot10^{19}}.
\]

The proof uses the upstream pair theorem and R1 finite certificates, the
Gaussian monotonicity and insertion identities already stated in
`GLOBAL_STRATEGY.md` and `QUANTITATIVE_GLOBAL.md`, and the new wide natural
neighborhood certificate. The latter proves that every normalized cone mass
vector with

\[
\|m-m_{\rm tree}\|_1\le\frac9{2500}=0.0036
\]

has total P4 cost at most `2*M_tree-10^-4`. This is the only new integration
input needed here. The constants below are checked by
`verify_stability_refinement.py`; the stored output is
`stability_refinement.json`. The wide inputs are
`certification/wide_neighborhood.json` and `certification/mass_neighborhood.json`.
The checker reads the exact upper cost endpoint and mass radius, checks the
required gap against the exact tree enclosure, and records the input hashes.

The resulting improvement remains far below the competitive record target.
The purpose of this note is to obtain a much stronger uniform consequence
from the existing rigidity method, while keeping its finite dependencies
and domain coverage explicit.

## 1. A sharper old-pair covariance estimate

Write `eta=17/(5*10^19)` and `v=sqrt(eta)`. As in the previous proof, if an old pair
cost is greater than `M_tree-eta`, the upstream R1 witnesses produce

\[
C_0\le C'\le C,\qquad
M_{\rm tree}-\Phi(C_0)\ge\lambda\|\theta-\theta_T\|_2^2,
\qquad\lambda=3349/10000.
\]

The nonlocal R1 alternatives are excluded because their certified deficit,
even after allowing the stated rounding allowance, exceeds `3*10^-12`.
The parameter lies in the upstream local box
`|theta-theta_T|_infinity<=10^-4`.

For a cross entry of the chart, write

\[
c(\theta)=(S/2+\sigma t_u)(S/2+\tau t_v),\qquad
S^2=K_0-a_{ab}(x+y),\qquad \sigma,\tau\in\{-1,1\}.
\]

On the whole straight parameter segment from `theta_T` to `theta`, one has
`3/2<S<2`, `a_ab<5`, and `|t_u|,|t_v|<=10^-4`. The first claim follows
directly from `S^2=4s-a_ab(x+y-2p)` and the same local-box bounds; it does not
require that every parameter on the segment be feasible. Thus

\[
|\partial_x c|,|\partial_y c|
\le\frac54\left(1+\frac1{7500}\right),\qquad
|\partial_{t_u}c|,|\partial_{t_v}c|\le1+10^{-4}.
\]

Consequently the squared gradient norm is at most

\[
L^2=2\left[\frac54\left(1+\frac1{7500}\right)\right]^2
       +2(1+10^{-4})^2<16\lambda.
\]

The two pair entries have gradient norm one and satisfy the same bound.
Cauchy-Schwarz along the parameter segment therefore gives

\[
\|C_0-C_T\|_{\max}<4\sqrt\eta.
\]

The analytic derivative estimate `partial_ij Phi<=-1/100` on the
radius-`10^-3` box around `C_T` is unchanged. Apply it to the initial
segment from `C_0` toward `C`: an increment `1/2000` would already lose
`1/200000` in cost, contradicting the assumed deficit. The full segment
therefore stays in that box, and

\[
\|C-C_0\|_{\max}\le100\eta.
\]

Since `100*sqrt(eta)<1`, the improved conclusion is

\[
\boxed{\|C-C_T\|_{\max}\le r:=5\sqrt\eta.}                 \tag{1}
\]

As in the original argument, a possibly singular endpoint is covered by
monotonicity after the positive-definite initial segment.

## 2. Row and star estimates with the useful constants retained

Assume both old pair costs exceed `M_tree-eta`, and put `u=sqrt(r)`.
The STAR and EC arguments from the previous proof give, at every core,

\[
0\le\bar k-p\le4r/3,\qquad
\sum_{e\text{ at the core}}(e-\rho)<7r,
\qquad |n_j/n_i-1|<5u
\]

on every required edge. Here `n_i=sqrt(G_ii)`, and the normalized masses
make `n_0=1`. Every label is at required-tree distance at most four from
label zero. For `u<=10^-4`,

\[
n_i\le(1+5u)^4<1003/1000,\qquad n_i^2<1007/1000.
\]

If a label has distance `d` from the root, its squared norm satisfies

\[
|G_{ii}-1|\le\frac{1002}{100}\,d\,u.                     \tag{2}
\]

For the upper endpoint, divide `((1+5u)^(2d)-1)` by `u`. This is a polynomial
with nonnegative coefficients, so it suffices to check `u=10^-4` for
`d=1,2,3,4`. The lower endpoint follows from Bernoulli's inequality.

The normalized star residual from the prior proof obeys

\[
\left\|\sum_{w\sim v}U_w-2\sqrt2 U_v\right\|
\le6u+8r.
\]

For the unnormalized rows `a_i=n_i U_i`, each adjacent norm difference is
at most `5*n_v*u`. Therefore

\[
R_v:=\sum_{w\sim v}a_w-2\sqrt2 a_v,
\qquad
\|R_v\|\le n_v(21u+8r)<22u.                             \tag{3}
\]

The last inequality is the exact scalar comparison
`(1003/1000)*(21+8*10^-4)<22`. This retains the information that the row
norms are already close to one, rather than bounding them by two.

## 3. Exact separation of every non-tree type

For a distance type `t`, define its raw-star energy

\[
E_t=\sum_{v=0}^3
\left(\sum_{w\sim v}g_t(w)-2\sqrt2 g_t(v)\right)^2.
\]

Every entry and every energy belongs to `Q(sqrt(2))`. Exact enumeration gives

\[
E_t=0\quad\text{for the six natural tree types},\qquad
E_t\ge E_0:=\frac{3-2\sqrt2}{16}\quad\text{for every other type}.        \tag{4}
\]

The minimum positive energy is attained at type indices
`4435,4437,4597,4599`. The verifier checks all 4,655 types using exact rational
coefficients and a rational enclosure of `sqrt(2)`; it checks equality
symbolically and strict inequalities with containing rational bounds.

Let `q_off` be the mass outside the six tree types. Source-coordinate
orthogonality, (3), and (4) imply

\[
E_0q_{\rm off}\le\sum_v\|R_v\|^2
<4\cdot22^2r.
\]

Since `4*22^2/E_0<200000`,

\[
\boxed{q_{\rm off}\le200000r.}                         \tag{5}
\]

This is a finite separation inequality over the entire cone. It does not
assume that the mass is supported near the tree face.

### 3.1 A sharper separation using normalized residuals

The preceding exact energy check is useful independently, but splitting the
off-face types gives the stronger bound

\[
\boxed{q_{\rm off}<2700r.}                              \tag{5a}
\]

First aggregate the STAR bounds over all four cores. For either old pair,
the sum of its two star means is at most
`(||Y||+C_ab+C_cd)/3 <= 2p+4r/3`. Thus

\[
\sum_v(\bar k_v-p)\le8r/3.
\]

Put `c0=sqrt((1+rho)/2)>0.98`. The derivative bound
`partial_i mu>=c0/4` now gives

\[
\sum_v\sum_{e\ni v}(e-\rho)<11r.
\]

At each core let `s_v=sum_i(c_i-c0)` and
`T_v=sum_i Y_i-(sum_i c_i)U_v`, in the notation of the previous proof.
Then

\[
\sum_v\|T_v\|^2\le16r,\qquad
\sum_v s_v<3r.
\]

The normalized star residual is `2*c0*T_v+K_v`, with `||K_v||<=4*s_v`.
Minkowski's inequality in the direct sum over cores therefore gives

\[
\sum_v\left\|\sum_{w\sim v}U_w-2\sqrt2U_v\right\|^2
\le(8\sqrt r+12r)^2<65r.                              \tag{5b}
\]

Now separate the off-face types into two classes.

**Types having a zero edge increment.** For a required edge `ij`, orient it
from a core label and write

\[
L_{ij}=\frac3{\sqrt2}G_{ij}-G_{ii}-G_{jj}.
\]

The per-type contribution vanishes when the increment is `+1` or `-1`.
When it is zero, it equals `(3/sqrt(2)-2)*g_t(i)^2`. On the other hand,

\[
L_{ij}\le\frac3{\sqrt2}n_i n_j(e_{ij}-\rho).
\]

Every core generator square is at least `1/8`, every norm product is less
than `1.007`, and the sum of required-edge excesses is bounded by the
core-incidence sum, hence is less than `11r`. Summing the nonnegative slack
over edges counts every zero-increment type at least once. Its total mass
`q_zero` satisfies

\[
q_{\rm zero}
\le\frac{8(3/\sqrt2)(1.007)(11r)}{3/\sqrt2-2}<1600r.
\]

The last comparison follows already from `3/sqrt(2)<15/7` and
`3/sqrt(2)-2>3/25`.

**Remaining non-tree types.** Their increments are all `+1` or `-1`, and
some core has a descending-neighbor count `d` other than one. At that core,

\[
\sum_{w\sim v}\frac{g_t(w)}{g_t(v)}-2\sqrt2
=\frac{d-1}{\sqrt2},
\]

whose absolute value exceeds `0.7`. In the normalized residual each neighbor
ratio is multiplied by `n_v/n_w`. This changes the displayed sum by at most
`5u*(3sqrt(2))<22u`. Since `g_t(v)^2>=1/8` and `n_v^2<1.007`, the squared
normalized residual contributed by this source coordinate, per unit mass,
is at least

\[
\frac{(0.7-22u)^2}{8(1.007)}>0.06
\]

for `u<=10^-4`. The source coordinates are orthogonal. By (5b), the total
mass of these remaining types is less than `65r/0.06<1100r`. Adding the two
classes proves (5a).

## 4. Recovering masses without the earlier norm overestimates

Use the same exact tree diagonal matrix `A` as before, on rows
`(0,1,2,3,4,8)`. These rows have depths `(0,1,2,3,1,4)`. The absolute column
sums of its inverse are

\[
(7/3,\ 5/3,\ 5/6,\ 5/12,\ 4/3,\ 1/6).
\]

Their depth-weighted sum is `79/12`. Equation (2) thus bounds the contribution
of the diagonal errors to the recovered tree masses by

\[
\frac{1002}{100}\frac{79}{12}u<66u.                     \tag{6}
\]

For any type `t`, let `d_t` denote its six selected diagonal entries. Exact
rational enumeration gives the much smaller combined off-face coefficient

\[
\max_{t\text{ off the tree face}}\left(1+\|A^{-1}d_t\|_1\right)=9/2.    \tag{7}
\]

To apply it, write `w` for the six retained tree masses and `w_T` for their
natural values. The exact diagonal identity is

\[
w-w_T=A^{-1}\left[(G_{ii}-1)_{i\in\{0,1,2,3,4,8\}}
                 -\sum_{t\text{ off}}m_td_t\right].
\]

The full L1 distance includes the off-face masses themselves. Combining
(5a)--(7), without renormalizing `w`, yields

\[
\boxed{\|m-m_{\rm tree}\|_1
\le66\sqrt r+\frac92q_{\rm off}
\le66\sqrt r+12150r.}                                 \tag{8}
\]

## 5. Closing the global dichotomy at `eta=17/(5*10^19)`

For the chosen `eta`, exact rational square comparisons give

\[
\sqrt\eta<\frac{729}{1250000000000},\quad
r<\frac{729}{250000000000},\quad
\sqrt r<\frac{27}{500000}.
\]

All preceding smallness conditions hold. Equation (8) gives

\[
\|m-m_{\rm tree}\|_1
<66\frac{27}{500000}+12150\frac{729}{250000000000}
=\frac{17997147}{5000000000}=0.0035994294<0.0036.
\]

The wide natural-neighborhood certificate therefore applies whenever both
old pair costs exceed `M_tree-eta`. It supplies total P4 deficit greater
than `10^-4`, which exceeds `eta`. If either old pair cost is at most
`M_tree-eta`, the nonnegative-gain decomposition immediately supplies deficit
at least `eta`. Thus

\[
M_4\le2M_{\rm tree}-\frac{17}{5\cdot10^{19}}
\]

holds on the complete closed simplex, including singular endpoints.
The P4-only structural reduction in `GLOBAL_REDUCTION.md` gives
`kappa<=M_4/4` uniformly over the required cubic kernels. Together with the
upstream circuit reduction, this yields

\[
\kappa\le\frac{M_{\rm tree}}2-\frac{17}{20\cdot10^{19}},
\qquad
c\ge1+\frac1{M_{\rm tree}-17/(10\cdot10^{19})}.
\]

The global deficit improves the previous explicit value by a factor
`3.4*10^35`. It remains a tiny improvement over the original method and does
not establish the much stronger record-scale P4-P7 target.
