# Natural-tree experiments for longer paths and short matching-score windows

## What this study establishes numerically

The earlier P4-P7 tree values are too large for the comparison target
`kappa=1341/10000=0.1341`. Larger paths have substantially smaller natural-tree
values. In this study P24 is already below that target, and P32-P128 have
larger numerical margins. A stationary seven-change window, involving eight
matching-anchor scores, is also below the target at the natural tree point.

These are **exploratory numerical expectations**, not rigorous integration
bounds or uniform inequalities over graph neighborhoods. They show that
longer paths and stationary matching-score windows are worth considering;
they do not establish a new cutwidth constant, a new circuit lower bound,
or an exact limiting value.

A subsequent [adversarial tree-face check](WINDOW_ADVERSARIAL.md) finds
that both the forward H7 bound and its reversal average fail away from
the natural point. The separate [window obstruction](WINDOW_OBSTRUCTION.md)
reduces that failure to a low-dimensional certificate. Thus the natural
window margins below are useful diagnostics, not evidence for a uniform
H7 upper bound over the entire projected local cone.

The final run uses eight independent scrambled Sobol sequences of
`2^20=1,048,576` Gaussian samples each, with seeds `88123,...,88130`.
All path lengths and windows share the samples within each scramble.
Two preliminary runs used disjoint seeds and were not included in the final
estimates. No optimized covariance or fitted control coefficient was chosen
using the final samples.

## 1. Longer-path results

The primary estimator uses the exact-mean sum of old P2 costs as a control
variate, explained below. Every Boolean cost is evaluated by a linear-time
dynamic program; no exponential cost table is created.

| Path vertices | Estimated natural-tree cost per vertex | Heuristic SE across scrambles | Raw direct mean | Raw heuristic SE |
|---:|---:|---:|---:|---:|
| 4 | 0.136133072 | 0.000008816 | 0.136104077 | 0.000035763 |
| 8 | 0.134890643 | 0.000007660 | 0.134856522 | 0.000021260 |
| 16 | 0.134264057 | 0.000004980 | 0.134260505 | 0.000018888 |
| 24 | 0.134055498 | 0.000003146 | 0.134048616 | 0.000013285 |
| 32 | 0.133955327 | 0.000001925 | 0.133960012 | 0.000011035 |
| 48 | 0.133853162 | 0.000002066 | 0.133846693 | 0.000008325 |
| 64 | 0.133803196 | 0.000002310 | 0.133804733 | 0.000007563 |
| 128 | 0.133724384 | 0.000000823 | 0.133727939 | 0.000005527 |

The P4 and P8 rows are checks against already studied small paths. The P4
reference has a separate interval certificate; these sampling estimates
are not used to alter that certificate.

The decrease through P128 suggests that endpoint costs are significant for
short paths. It does not prove convergence at a stated rate, establish an
exact limiting constant, or show that the natural point maximizes a larger
path cone. A uniform bound over an entire relaxed path domain must in
particular handle its natural-tree point, but passing this necessary check
is far from sufficient.

## 2. Stationary matching-anchor windows

Let `b_i` be the sign of the **single pendant matching-anchor score** at
backbone vertex `i` in the infinite natural tree process. These are not
the two-endpoint boundary arrays used for finite Pk costs. Their stationary
Gaussian covariance is

\[
K_{ii}=1,\qquad
K_{ij}=\left(p_0+\frac{|i-j|}{4}\right)2^{-|i-j|/2}\quad(i\ne j),
\qquad p_0=\frac{1+2\sqrt2}{4}.
\]

Set `D_i=1[b_i != b_(i+1)]`. For odd `p`, the sampled window function is

\[
H_p(i)=\sum_{j=0}^{p-1}(-1)^j\prod_{a=0}^{j}D_{i-a}.
\]

It takes values in `{0,1}`. Expanding the greedy matching recurrence over
a window shows why the odd truncation is useful: it is the value obtained
with incoming state zero and bounds the final state for either possible
incoming state. Turning this local identity into a universal cutwidth
certificate requires its separate graph/ordering argument and a uniform
Gaussian bound. This numerical study evaluates only its natural-tree
expectation.

For each Gaussian sample, the script averages over every valid window in
the 128-core pendant sequence: `128-p` positions. The positions within a
sample are correlated; they are **not** treated as independent samples.

| Changes in window | Anchor-score coordinates | Estimated expectation with P2 control | Heuristic SE | Raw direct mean | Raw heuristic SE |
|---:|---:|---:|---:|---:|---:|
| 3 | 4 | 0.137694477 | 0.000001935 | 0.137698031 | 0.000006523 |
| 5 | 6 | 0.134044743 | 0.000001419 | 0.134048297 | 0.000005848 |
| 7 | 8 | 0.133686826 | 0.000001668 | 0.133690381 | 0.000005076 |
| 9 | 10 | 0.133652458 | 0.000002163 | 0.133656013 | 0.000005476 |

The H7 estimate is approximately `0.00041317` below the comparison target
`0.1341`. H5 has a smaller numerical margin, approximately `0.00005526`.
H9 is lower again, at the cost of a larger local Gaussian functional.
These margins concern only the sampled natural point.

The JSON also stores an alternative sign-change control. Its means for
H3, H5, H7, and H9 are respectively approximately `0.137697092`,
`0.134045821`, `0.133686094`, and `0.133648901`, with heuristic standard
errors between `3.1e-6` and `4.1e-6`. Both controls have preselected coefficient
one; neither was fitted to the observed samples.

## 3. Common Gaussian process and Boolean evaluation

The full field contains 257 edge scores: the left entering backbone edge,
and for every core vertex `v=0,...,127`, its pendant edge and forward
backbone edge. The common covariance is built from the natural tree vertex
correlation

\[
r_d=(1+d/3)2^{-d/2}
\]

by summing the four endpoint correlations and dividing by `2+2*r_1`.
This is the same construction used in the earlier common-tree study.

In the stored coordinate order, the entering edge is index zero, pendant
edge `v` is index `2v+1`, and forward edge `v` is index `2v+2`. A Pk boundary
therefore uses

```
[0] + [2*v+1 for v in range(k)] + [2*k].
```

The code checks every requested marginal against the separately constructed
`tree_covariance(k)`. The maximum observed discrepancy was less than
`1e-12` (the recorded per-length values are approximately machine precision).

A floating-point eigendecomposition gives the expected rank 130. Modes at
most `1e-11` are removed, and the reconstructed common covariance differs
from its input by about `4.22e-15` entrywise. These are floating-point
consistency checks, not interval spectral certificates. Normal coordinates
come from inverse-normal transforms of the scrambled Sobol points.

For each requested path, `path_gadgets.cost` computes the exact integer
minimum split-vertex count in `O(k)` time per sample. Thus the experiment
does not rely on the exponentially sized Boolean table limited to small
paths in earlier scripts.

## 4. Exact-mean control variates

For even `k`, let `B_k` be the sum of the optimized P2 costs on consecutive
vertex pairs `(0,1),(2,3),...`, keeping their connecting backbone-edge
scores fixed at their sampled values. Each old pair has exactly the natural
P2 marginal, so

\[
\mathbb E B_k=\frac{k}{2}M_{\rm tree}.
\]

For a pair with boundary bits `a,b,c,d`, its Boolean cost is evaluated as

\[
\mathbf1[a\ne b]+\mathbf1[c\ne d]
+\mathbf1[a=b,\ c=d,\ a\ne c].
\]

Joint path optimization can only improve this feasible separate-pair
assignment. The program checks `B_k-F_k>=0` for every sample and uses

\[
\widehat{\kappa}_k
=\frac{M_{\rm tree}}2
-\frac1{Nk}\sum_{s=1}^N(B_k^{(s)}-F_k^{(s)}).
\]

The exact tree enclosure supplies the control mean; the tiny width of that
enclosure is negligible on the sampling scale but is retained in the JSON.

For a window, write `bar H_p` for its average over valid positions within
one 128-core sample. Stationarity makes its expectation equal to that of
any one window. The primary estimator is

\[
\widehat h_p^{\rm pair}
=\frac{M_{\rm tree}}2
+\frac1N\sum_{s=1}^N
\left(\overline H_p^{(s)}-\frac{B_{128}^{(s)}}{128}\right).
\]

An additional control uses the adjacent matching-score correlation
`s=(2+sqrt(2))/4`, hence

\[
\mathbb E D_i=\frac{\arccos s}{\pi}.
\]

Subtract the average of `D_i` over the same window-ending positions, then
add this known expectation. Both control identities preserve the desired
expectation in exact arithmetic. Gaussian generation, factorization, and
transcendental evaluations in this experiment use ordinary floating point.

## 5. Uncertainty, reproduction, and scope

For each reported estimator, the standard error is the sample standard
deviation of its **eight replicate estimates**, divided by `sqrt(8)`.
It reflects variation between independent scrambles. It is not a guaranteed
confidence interval, a deterministic quadrature remainder, or a proof of
the ordering of nearby expectations. Raw direct estimates and their
replicate variation are included to make the variance reduction transparent.

The final run processed 8,388,608 common Gaussian samples in batches of
8,192. It took approximately 61 seconds in the development environment.
Reproduce it from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/long_path_study.py \
  --sizes 4 8 16 24 32 48 64 128 \
  --power 20 --repeats 8 --seed 88123 --batch-power 13 \
  --output research/gadget_insertion/long_path_study.json
```

Requirements are Python, NumPy, and SciPy. The JSON records software
versions, seeds, sample counts, batching, all replicate path values, raw
means, both window control estimates, covariance diagnostics, and the exact
enclosure used for the P2 control mean.

Natural-tree behavior leaves room below `0.1341` for sufficiently long
paths and for H5 or H7 matching-score windows. The subsequent adversarial
study shows that this headroom does not extend to a bare uniform H7 bound
over the projected tree face, even after reversal averaging. A competitive
theorem therefore needs further structure: for example, a bound using
consistency between overlapping windows, a cycle-wide average, or an
alternative rule. A uniform longer-path bound remains a separate question.
No such theorem is asserted by these experiments.
