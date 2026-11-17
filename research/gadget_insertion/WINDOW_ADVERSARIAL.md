# Adversarial checks for the seven-change matching-score window

## Main finding

The natural-tree expectation of the seven-change window is numerically
below `0.1341`, as reported in [LONG_PATH_STUDY.md](LONG_PATH_STUDY.md).
That headroom does **not** persist throughout the projected P8 tree face.
A limited search finds feasible positive tree-face masses with estimated
forward-window expectation about `0.143893`, and other masses with
reversal-averaged expectation about `0.145481`. A candidate at full mass
L1 distance about `0.0112551` from the natural point already exceeds the
target for both objectives.

The numerical results here do not certify integration errors. They identify
an obstruction to a proposed uniform local bound. The separate
[rigorous certificate](WINDOW_OBSTRUCTION.md) uses a positive rational
mass vector and proves that both orientations and their average have
expectation strictly above `17/125=0.136`, without integrating an
eight-dimensional sign functional. None of these relaxed mass points is
asserted to arise as a particular window in a finite cubic graph. Thus the
finding concerns uniform bounds over this relaxation, and does not rule
out a bound that uses consistency between overlapping windows or averages
over an entire cycle.

## 1. Precisely which function and domain were searched

There are eight path core vertices, with one external matching-anchor
score selected for each core vertex. The ten boundary scores of P8 have
owners

```
[0, 0, 1, 2, 3, 4, 5, 6, 7, 7].
```

The projection keeps boundary coordinates `[1,2,3,4,5,6,7,8]`, whose owners
are `[0,1,2,3,4,5,6,7]`. The covariance is constructed by the existing
`path_cone.covariance` function, including both vertex and edge-sum
normalizations, and then restricted to those eight coordinates.

Let the resulting signs be `b_0,...,b_7`, and put
`D_i=1[b_i != b_(i+1)]` for `i=0,...,6`. The forward objective is

\[
H_7^F=D_6-D_6D_5+D_6D_5D_4-\cdots+D_6D_5\cdots D_0.
\]

The reverse objective is the same function after replacing `b_i` by
`b_(7-i)`. The second searched objective is
`H_7^sym=(H_7^F+H_7^R)/2`. An individual window can have different forward
and reverse expectations away from the natural stationary point; the
averaging is therefore performed explicitly. The exact tables have 256
entries, with values in `{0,1}` for the forward function and `{0,1/2,1}`
for its symmetrization. Both functions are invariant under complementing
all eight signs.

Only the ten types returned by `path_cone.tree_support(8)` are used. Their
ordered rows are stored in the JSON, so there is no reliance on a full
enumeration of the P8 cone. The natural masses in that order are

\[
m_T=\left(\frac13,\frac13,\frac16,\frac1{12},\frac1{24},
\frac1{48},\frac1{96},\frac1{192},\frac1{384},\frac1{384}\right).
\]

The numerical search parametrizes positive masses by

\[
m_j(z)=\frac{(m_T)_j e^{z_j}}
 {\sum_{a=0}^{9}(m_T)_a e^{z_a}},
\qquad z_9=0,\quad -0.8\le z_j\le0.8\quad(0\le j<9).
\]

This is a bounded subset of the tree-face simplex. It does not contain all
positive masses on that face, much less the whole P8 cone.

## 2. Numerical results

Each row below is evaluated on six holdout scrambles, each containing
`2^18=262,144` residual Gaussian samples. The same holdout points are used
for every candidate. The heuristic standard error is the sample standard
deviation of the six replicate estimates, divided by `sqrt(6)`.

| Candidate | Full mass L1 distance from natural | Forward H7 mean | Forward heuristic SE | Symmetrized H7 mean | Symmetrized heuristic SE |
|---|---:|---:|---:|---:|---:|
| Natural | 0 | 0.133685048 | 0.000015511 | 0.133684801 | 0.000007624 |
| Forward fit, natural start | 0.491377574 | 0.143810218 | 0.000015172 | 0.126971451 | 0.000008644 |
| Forward fit, best probe start | 0.011255086 | 0.143892840 | 0.000015318 | 0.138804048 | 0.000008800 |
| Symmetric fit, natural start | 0.664137755 | 0.143598997 | 0.000018057 | 0.145480629 | 0.000008601 |
| Symmetric fit, best probe start | 0.664140155 | 0.143598998 | 0.000018057 | 0.145480628 | 0.000008601 |

The natural row agrees on the numerical scale with the independently
sampled long-path study. The candidate called `forward_best_probe` has
paired estimated increases over the natural row of about `0.010207792`
forward and `0.005119247` symmetrized. Their paired heuristic standard
errors are about `0.000016182` and `0.000008579`, respectively. The common
holdout points are used when computing these paired differences.

Its masses, in the stored tree-face order, are approximately

```
[0.3352289267584118,    0.3352289267584118,
 0.1676144633792059,    0.08380723168960295,
 0.04190361584480148,   0.02095180792240074,
 0.01047590396120037,   0.0010575242663063858,
 0.002554817650889257,  0.0011767817687694149].
```

The JSON preserves all candidate masses, covariances, training values,
holdout replicates, paired differences, optimizer status, and floating
minimum eigenvalues. The latter are positive for every reported
candidate, but are consistency diagnostics rather than interval
certificates.

## 3. Integration and search procedure

For an eight-dimensional positive definite covariance `C`, set

\[
\alpha=\frac{0.95}{\mathbf1^TC^{-1}\mathbf1},
\qquad R=C-\alpha\mathbf1\mathbf1^T.
\]

In exact arithmetic `R` is positive definite. Write the Gaussian vector
as `X=sqrt(alpha)*Z*1+Y`, with independent standard normal `Z` and
`Y~N(0,R)`. For a fixed residual sample `Y`, the eight thresholds
`-Y_i/sqrt(alpha)` divide the scalar `Z` axis into intervals on which the
sign pattern is constant. Summing the exact table value times the normal
probability of each interval integrates `Z` analytically. The table values
outside the first and last thresholds are zero because all signs agree.

The remaining integral is sampled using an eight-dimensional scrambled
Sobol sequence, inverse-normal transformation, and the Cholesky factor of
`R`. This representation varies smoothly with positive definite `C`
apart from harmless threshold-order changes. Matrix operations and normal
CDF evaluations use ordinary floating point. In particular, this
conditional integration is a variance-reduction method, not a certified
quadrature rule.

The search uses `2^13=8192` training points with Sobol seed `99117`. There
are 83 initial probes: the natural point, 64 clipped Gaussian log-mass
perturbations using random seed `99119`, and 18 single-coordinate
perturbations of size `0.35`. For each objective, L-BFGS-B is run from
both the natural point and the best probe for that objective. Each fit
uses the fixed box above, at most 80 iterations, and fixed finite
difference step `1e-5`. All four stored fits report numerical convergence;
this does not certify a global or local maximum of the exact objective.

The holdout seeds are `99231,...,99236`, disjoint from the training seed
and from the long-path experiment. They give 1,572,864 residual Gaussian
samples per candidate. All four fits are reported, rather than retaining
only the fit with the largest holdout value. Variation across scrambles
is reported as a heuristic standard error, not a rigorous confidence
bound or deterministic integration remainder.

## 4. A low-dimensional route to a rigorous obstruction

The Boolean identity

\[
H_7^F=D_6(1-D_5)
 +D_6D_5D_4(1-D_3)
 +D_6D_5D_4D_3D_2(1-D_1)
 +D_6D_5D_4D_3D_2D_1D_0
\]

implies the pointwise bound `H_7^F >= D_6(1-D_5)`. If `C` is the
correlation matrix of the eight selected Gaussian scores, the elementary
Gaussian sign identity gives

\[
\mathbb E[D_6(1-D_5)]
=\frac14+\frac{\arcsin C_{56}-\arcsin C_{67}-\arcsin C_{57}}{2\pi}.
\]

The reverse lower bound is

\[
\mathbb E[D_0(1-D_1)]
=\frac14+\frac{\arcsin C_{21}-\arcsin C_{10}-\arcsin C_{20}}{2\pi}.
\]

Averaging these bounds gives a lower bound for the symmetrized objective.
At the numerical `symmetric_natural` candidate, the two lower bounds are
approximately `0.1326800743` and `0.1359750832`; their average is about
`0.1343275787`, already above `0.1341`. This reduces the obstruction
certificate to interval evaluation of a handful of covariance entries
and arcsines, with no eight-dimensional Gaussian integration.

The resulting separate certificate is documented in
[WINDOW_OBSTRUCTION.md](WINDOW_OBSTRUCTION.md). It uses a simpler positive
rational mass vector on the same tree face:

\[
m=\frac1{1140}(128,640,64,160,80,40,20,2,5,1).
\]

It certifies the forward three-score probability in
`[0.136076509488866,0.136076509488875]`, the reverse probability in
`[0.136076509488864,0.136076509488876]`, and their average in
`[0.136076509488865,0.136076509488875]`. These are certified lower bounds
for the corresponding full H7 expectations, not enclosures of the full
H7 expectations themselves. Every lower endpoint exceeds `0.136`.
This witness is outside the bounded log-mass search box used for the
exploratory optimization, while remaining inside the same positive
tree-face simplex.

## 5. Reproduction and implication for further work

From the repository root, run

```sh
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/window_adversarial.py \
  --output research/gadget_insertion/window_adversarial.json
```

The stored run took approximately 11 seconds. It requires Python, NumPy,
and SciPy. Its source, tables, seeds, exact natural masses, candidate
covariances, and replicate estimates are included in the repository.

A proof based on a uniform `E H7 <= 0.1341` over this whole projected
tree face cannot follow from natural-point behavior. Reversal averaging
alone does not remove the obstruction. A competitive window proof would
need additional structure, such as an inequality for the sum around a
cycle, a telescoping potential that absorbs asymmetric local costs, or
constraints enforcing which local distributions can occur together. The
natural long-path estimates remain useful for assessing those directions,
but they do not supply the missing uniform bound.
