# A larger mass neighborhood from interval derivatives

## Certified result

Let `m_T` be the exact natural P4 mass vector in the full 4,655-type
simplex. The finite interval calculation in `mass_neighborhood.py` proves
that every normalized nonnegative mass vector satisfying

\[
\|m-m_T\|_1\le\frac9{2500}=0.0036
\]

has boundary correlation matrix satisfying

\[
C_4(m)_{ij}\ge C_4(m_T)_{ij}-\frac1{300}\qquad(i<j).
\]

The largest certified downward entry change is the exact rational

\[
\frac{18372519041787927}{5629499534213120000}
<\frac1{300}.
\]

Therefore the independently integrated comparison point from
`WIDE_NEIGHBORHOOD.md` applies throughout this full-simplex ball:

\[
\Phi_{P_4}(C_4(m))<2M_{\rm tree}-\frac1{10000}.
\]

This result allows arbitrary off-support masses. It is a local theorem and
does not assert that the whole cone satisfies this cost bound. The new
calculation performs no Gaussian quadrature; it checks derivatives of the
explicit two-normalization map, then invokes the existing wide comparison
certificate.

## 1. Signed mass changes and the containing Gram box

For each distance type `t`, put

\[
K_t(i,j)=g_t(i)g_t(j)=2^{-(t_i+t_j)/2},\qquad
G(m)=\sum_t m_tK_t.
\]

Every coefficient belongs to `Q(sqrt(2))`. The implementation constructs
containing intervals from exact dyadic constants and an outward square root
of two. It checks all 4,655 type coordinates and verifies the exact rational
identity `G(m_T)_ii=1` at each of the ten vertex labels.

Write `h=m-m_T`, and let `S` be the six-element support of `m_T`. Since
`sum h=0`, its positive and negative total masses both equal `||h||_1/2`.
Moreover, every negative coordinate of `h` lies in `S`, because `m_T` is
zero elsewhere and `m` is nonnegative. Consequently, for any scalar
coefficient array `a_t`,

\[
\frac{\|h\|_1}{2}
 (\min_t a_t-\max_{t\in S}a_t)
\le\sum_t h_ta_t\le
\frac{\|h\|_1}{2}
 (\max_t a_t-\min_{t\in S}a_t).
\]

Apply this entrywise to `K_t`. Since each lower range is nonpositive and
each upper range is nonnegative, replacing `||h||_1` by the chosen radius
gives a containing Gram box. Every point of the segment
`m_T+s h`, `0<=s<=1`, is in the simplex and in the same ball, so its Gram
matrix lies in that box. The program verifies positive lower bounds for
all Gram diagonals.

## 2. Derivatives through both normalizations

The normalized vertex covariance is

\[
V_{ab}=\frac{G_{ab}}{\sqrt{G_{aa}G_{bb}}}.
\]

For a fixed matrix direction `D`, its derivative is exactly

\[
\dot V_{ab}=\frac{D_{ab}}{\sqrt{G_{aa}G_{bb}}}
-\frac{V_{ab}}2
 \left(\frac{D_{aa}}{G_{aa}}+\frac{D_{bb}}{G_{bb}}\right).
\]

For `a=b`, this derivative is identically zero, which the code preserves
exactly. For the boundary edges `e_i=(o_i,4+i)`, define

\[
W_{ij}=\sum_{a\in e_i}\sum_{b\in e_j}V_{ab},\qquad
C_{ij}=\frac{W_{ij}}{\sqrt{W_{ii}W_{jj}}}.
\]

Here `W_ii=2+2V_(o_i,4+i)` and
`dot W_ii=2 dot V_(o_i,4+i)`. Applying the same normalization formula a
second time gives `dot C_ij`. The program encloses these expressions over
the entire Gram box for each of the 4,655 directions

\[
D_t=K_t-G(m_T).
\]

Subtracting the common center direction is valid: since `sum h=0`,
`sum h_t D_t=sum h_t K_t`. Thus it does not change the derivative along
the actual mass displacement. This shift is used only to reduce interval
width.

Every actual vertex and boundary covariance has positive entries and is
positive semidefinite. Intersecting their interval value enclosures with
`[0,1]` is therefore valid even if unrelated matrices in the containing
entry box would fail to be positive semidefinite. No derivative is clipped;
the intersections apply only to values with the stated mathematical bound.

## 3. A one-sided derivative bound along the entire mass segment

Let `a_t(s)` be the derivative of one selected boundary correlation at
`m_T+s h` in direction `D_t`. The interval computation bounds every such
number for all `0<=s<=1`. The signed-mass observation from Section 1 gives

\[
\frac{d}{ds}C_{ij}(m_T+s h)
\ge-\frac{\|h\|_1}{2}
 \left(\max_{t\in S}a_t(s)-\min_t a_t(s)\right).
\]

For each of the fifteen correlations, the script takes the upper interval
endpoint of half this difference, then multiplies by the exact rational
radius. Integrating the derivative over the unit segment proves the
advertised downward-change bound. The maximum over all fifteen entries
is strictly smaller than `1/300` by exact rational comparison.

The region is deliberately one-sided: a covariance entry increasing above
the natural value does not harm the subsequent submodular Gaussian
comparison. Bounding an absolute change by an unrelated worst direction
would lose this information.

## 4. Reproduction and dependencies

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -B \
  research/gadget_insertion/certification/mass_neighborhood.py \
  --output research/gadget_insertion/certification/mass_neighborhood.json
```

The default radius is `9/2500`. The run takes less than a second in the
development environment. Its per-pair slope bounds and exact rational
consequences are saved in the JSON file.

The arithmetic trust base is the pinned upstream `ivarith.py`, Git blob
`b99d01f757f2b3d75151ba9403a547c6b8c60477`. The script checks that hash.
Dyadic constants are produced by `ldexp` with small integer exponents;
all subsequent arithmetic and square roots use containing intervals.
The assumptions are the same documented IEEE binary64 assumptions as in
the previous interval certificates.

The separate comparison-point input is `wide_neighborhood.json`, whose
Gaussian integral can be replayed with `wide_neighborhood.py`. The mass
checker verifies its stated comparison radius and compares its exact cost
upper endpoint against the independent exact lower enclosure of
`M_tree`. Input hashes identify the source and certificate versions used.
