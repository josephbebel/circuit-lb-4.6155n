# A rigorous obstruction to the relaxed H7 window bound

## Result and scope

There is an explicit positive rational mass vector on the ten tree types
of the P8 distance-type relaxation for which

\[
\mathbb E H_7^{\rm forward}>\frac{17}{125}=0.136,
\qquad
\mathbb E H_7^{\rm reverse}>\frac{17}{125},
\qquad
\mathbb E\frac{H_7^{\rm forward}+H_7^{\rm reverse}}2>
\frac{17}{125}.
\]

Thus the pointwise upper target `1341/10000=0.1341` fails on this relaxed
face for either orientation and for their average. This is a certified
obstruction, obtained from three-score marginal probabilities. No
high-dimensional integration or sampling is used in the certificate.

**This is a statement about the relaxed covariance face. No cubic graph
realizing the rational mass vector is supplied.** It does not disprove a
window bound under additional graph-realizability or compatibility
constraints, a bound averaged over overlapping windows, or a different
long-path insertion objective. It proves no lower bound on cubic-graph
cutwidth and no new circuit lower bound.

The experimental search that motivated the certificate is documented in
`WINDOW_ADVERSARIAL.md`. The simple rational witness below is a separate
point on the same ten-type face; it need not lie within that experiment's
restricted relative log-mass search box.

## 1. An exact three-score lower bound

For eight bits `b_0,...,b_7`, write

\[
D_i=\mathbf1[b_i\ne b_{i+1}],\qquad
H_7^{\rm forward}=
\sum_{\ell=1}^{7}(-1)^{\ell+1}
\prod_{j=0}^{\ell-1}D_{6-j}.
\]

The nested products imply `H_7` is either zero or one. If `D_6=1` and
`D_5=0`, precisely its first product is nonzero, so `H_7=1`. Therefore

\[
H_7^{\rm forward}\ge A_f:=D_6(1-D_5),\qquad
H_7^{\rm reverse}\ge A_r:=D_0(1-D_1).                 \tag{1}
\]

Adding these inequalities also bounds their orientation average. The
verifier independently checks (1), the average inequality, and the Boolean
range on all 256 assignments.

Let `X` be a centered Gaussian vector with unit coordinate variances and
correlation matrix `C`, and set `b_i=1[X_i<0]`. Write `s_i` for its sign
in `{−1,1}`; zero scores have probability zero. Then

\[
A_f=\frac{(1+s_5s_6)(1-s_6s_7)}4
=\frac{1+s_5s_6-s_6s_7-s_5s_7}{4}.
\]

The Gaussian pair-sign identity
`E[s_i s_j]=(2/pi) asin(C_ij)` follows from the rotationally invariant
two-dimensional Gaussian: the sign-disagreement sectors have total angle
`2 acos(C_ij)`. This is also Lemma G(c) of the repository's `proof.md`.
Consequently

\[
q_f:=\mathbb E A_f
=\frac14+
\frac{\arcsin C_{56}-\arcsin C_{67}-\arcsin C_{57}}{2\pi},
\]

\[
q_r:=\mathbb E A_r
=\frac14+
\frac{\arcsin C_{12}-\arcsin C_{01}-\arcsin C_{02}}{2\pi}.             \tag{2}
\]

Equations (1) and (2) reduce the required lower certificate to six
correlation entries, square roots, and arcsines.

## 2. The exact rational covariance witness

Use the extended P8 tree with core vertices `0,...,7` and ten boundary
leaves whose owners, in order, are

```
[0, 0, 1, 2, 3, 4, 5, 6, 7, 7].
```

Leaf `j` has vertex label `8+j`. Let `d` denote ordinary tree distance and
define the corresponding eighteen-coordinate ray by

\[
t^{(j)}_v=d(v,8+j)-d(0,8+j),\qquad
u^{(j)}_v=2^{-t^{(j)}_v/2}.
\]

These are exactly the ten branch-ordered rays returned by
`path_cone.tree_support(8)`. The certificate constructs them independently
from integer tree distances. It checks the edge distance constraints and
the generic descent condition at every core vertex, without enumerating
the full P8 cone.

In this branch order, take

\[
m=\frac1{1140}
(128,640,64,160,80,40,20,2,5,1).                       \tag{3}
\]

All masses are positive and their sum is exactly one. For comparison, the
natural tree masses have numerators
`(128,128,64,32,16,8,4,2,1,1)` over `384`. The witness reweights those
branches by factors `(1,5,1,5,5,5,5,1,5,1)` and normalizes.

Form the positive semidefinite Gram matrix and its vertex correlation
matrix:

\[
G=\sum_{j=0}^{9}m_j u^{(j)}(u^{(j)})^T,
\qquad
V_{xy}=\frac{G_{xy}}{\sqrt{G_{xx}G_{yy}}}.
\]

The matching-anchor boundary indices are `1,...,8`, with one anchor at
each core vertex. If index `i` has endpoints `(a_i,c_i)`, form

\[
E_{ij}=V_{a_i a_j}+V_{a_i c_j}+V_{c_i a_j}+V_{c_i c_j},
\qquad
C_{ij}=\frac{E_{ij}}{\sqrt{E_{ii}E_{jj}}}.
\]

This is the same two-stage normalization as the original distance-type
relaxation. The diagonal entries of `G` are positive. The projected
eight-score matrix `C` is furthermore certified positive definite by
interval Cholesky with strictly positive pivots.

## 3. Certified numerical enclosures

The verifier produces the following outward decimal intervals. Each JSON
record also retains the exact rational values of the binary64 endpoints.

| Quantity | Certified enclosure |
|---|---|
| `q_f` | `[0.136076509488866, 0.136076509488875]` |
| `q_r` | `[0.136076509488864, 0.136076509488876]` |
| `(q_f+q_r)/2` | `[0.136076509488865, 0.136076509488875]` |

Every lower endpoint strictly exceeds `17/125`, which in turn strictly
exceeds `1341/10000`. By (1), these are lower bounds for the respective
H7 expectations. They are not claimed to be the full H7 expectations or
global maxima over the relaxed face.

The arithmetic source is the original `mtree/enclosure/ivarith.py`, pinned
before import to SHA-256

```
cb1aec517a74efc071374e5cd4426e9bd05958e8d03559ed1d657d6e85186145
```

Its elementary operations round outward. Its arcsine endpoints are checked
using `mpmath` interval sine. The verifier independently checks the inherited
pi enclosure against `mpmath.iv.pi`, converting its binary64 endpoints
through exact dyadic rationals. Every rational input is first enclosed by
exact comparisons with the binary64 endpoint fractions. The covariance
normalizations, Cholesky pivots, probability formulas, and final comparisons
all use intervals. An arithmetic-source hash mismatch aborts verification.

## 4. Replay and consequence for further work

From the repository root:

```sh
python -B research/gadget_insertion/certification/window_obstruction.py \
  --output research/gadget_insertion/certification/window_obstruction.json
```

The run takes less than one second in the development environment and
requires Python, NumPy, and mpmath. It records the exact rational support,
all projected covariance intervals, Cholesky pivots, Boolean verification,
probability enclosures, source hash, and runtime versions.

The finite-window conditional reduction in `CYCLE_WINDOWS.md` remains
valid. What this witness rules out is proving its competitive pointwise
H7 hypothesis merely over the entire projected P8 tree face. Averaging the
two window orientations does not fix that relaxed-domain obstruction.
Further progress must use a restriction justified by actual graph fields,
constraints coupling multiple windows, or an improved objective.
