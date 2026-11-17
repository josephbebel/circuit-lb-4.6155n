# A finite-window objective for joint path and cycle insertion

## Purpose and scope

The refreshed frontier obstructs the current uniform P4 objective; see
`FRONTIER_OBSTRUCTION.md`. Joint optimization along longer paths or whole
2-factor cycles changes that objective. This note gives an exact Boolean
reduction and rigorous finite-window bounds for the changed rule.

The Boolean statements are proved below and checked by `cycle_windows.py`.
Section 6 gives a conditional bounded-gadget reduction with explicit endpoint
overhead, and Section 7 gives a compatible telescoping-potential interface.
Section 8 records a proved obstruction to the bare H7 target on the relaxed
tree face. No uniform Gaussian upper bound or new circuit coefficient is
asserted by this note.

## 1. Split vertices form a vertex cover of the change edges

Consider a path or cycle whose edges are free to receive side bits. Each
vertex `i` has one fixed external side bit `b_i`. A vertex has cost zero
precisely when every incident free-edge bit equals `b_i`; otherwise it has
cost one. For each path or cycle edge, put

\[
D_i=\mathbf1[b_i\ne b_{i+1}].
\]

If `D_i=1`, at least one endpoint must be split: the shared free edge cannot
agree with both boundary bits. Thus the split vertices form a vertex cover
of the graph consisting of change edges.

Conversely, take any vertex cover `S` of the change edges. Give every free
edge incident to a vertex outside `S` that vertex's boundary bit. This is
consistent if both endpoints lie outside `S`, since their bits must agree.
Assign remaining free edges arbitrarily. Every vertex outside `S` is
unsplit, so the cost is at most `|S|`. Hence the optimum equals the minimum
vertex-cover size of the change-edge graph.

For a path the change-edge graph is a disjoint union of paths. For a cycle
it is either a union of paths or the whole cycle. The latter case is possible
only for an even cycle, because the boundary bits alternate around it.
Its minimum vertex-cover size equals its maximum matching size in every
case. This also follows directly: a path with `r` edges has minimum cover
and maximum matching size `ceil(r/2)`, and an even cycle has both equal to
half its order.

For the actual cubic path `P_k`, the two endpoints each have two fixed
boundary bits. If an endpoint's bits disagree, it is forced to split and
contributes one. Remove such forced endpoints from the change-edge graph.
At an endpoint whose bits agree, use their common value as its color.
The exact optimum is

\[
\#\{\text{forced endpoints}\}
+\nu(\text{remaining change-edge graph}).
\]

Thus the difference from the one-boundary-bit path mechanism is confined to
the two endpoints. This identity is checked against the existing exact cubic
path cost tables.

## 2. Greedy matching and the alternating expansion

On a path, scan the change edges in order, greedily selecting a present edge
unless the preceding edge was selected. Write its selection indicator as
`z_i`. Then

\[
z_i=D_i(1-z_{i-1}),\qquad z_{-1}=0.
\]

This chooses a maximum matching on each component of consecutive change
edges. Expanding the recurrence through `K` terms gives

\[
z_i=\sum_{\ell=1}^K(-1)^{\ell+1}
       \prod_{j=0}^{\ell-1}D_{i-j}
 +(-1)^K\left(\prod_{j=0}^{K-1}D_{i-j}\right)z_{i-K},
\]

with nonexistent path edges taken to have `D=0`. Therefore odd truncations
are upper bounds and even truncations are lower bounds, pointwise.

For a cycle with a missing change edge, start after such an edge and use the
same recurrence around the cycle. It closes consistently. For a cycle
whose change edges are all present, the cycle is even and either alternating
perfect matching supplies a cyclic solution of the recurrence. Thus the
same pointwise truncation bounds apply on cycles, including wraparound.

## 3. A fixed finite-window upper bound

For a window of `K+1` boundary bits, define

\[
H_K(b_0,\ldots,b_K)
=\sum_{\ell=1}^K(-1)^{\ell+1}
 \mathbf1[b_{K-\ell},\ldots,b_K\text{ alternate}].
\]

Nested alternating-block indicators make `H_K` a Boolean-valued observable:
its value is always zero or one. It is invariant under complementing all
bits. For odd `K`, the exact cycle split optimum satisfies

\[
F_{\rm cycle}(b)
\le\sum_{i\in V(\rm cycle)}H_K(b_{i-K},\ldots,b_i),
\]

where indices are cyclic. No probabilistic independence or covariance
assumption enters this inequality. Moreover,

\[
0\le\sum_i H_K-F_{\rm cycle}
\le\sum_i\mathbf1[b_{i-K},\ldots,b_i\text{ alternate}].
\]

Taking expectations under any random boundary scores preserves these
inequalities. For a stationary process, put

\[
q_\ell=\Pr(b_0,b_1,\ldots,b_\ell\text{ alternate}).
\]

The expected odd-window cost per vertex is exactly

\[
\mathbb E H_K=\sum_{\ell=1}^K(-1)^{\ell+1}q_\ell.
\]

The infinite alternating series is useful only after its limiting and
boundary terms are justified. The finite odd-window inequality itself
already gives a valid upper objective and requires no limiting argument.

## 4. Exact tables and a necessary warning about submodularity

`cycle_windows.py` constructs the tables for `K=1,3,5,7,9`, their exact
integer Walsh numerators, their maximum Walsh degrees, and every mixed
Boolean difference. Complement symmetry makes all odd-degree coefficients
zero. The six-variable table `H_5` and eight-variable table `H_7` offer
concrete finite observables for a Gaussian calculation without enumerating
an entire long-path sign table.

For `K>=3`, the individual window observable is generally **not
submodular**. It therefore does not automatically inherit the monotone
Gaussian comparison used by the P4 cost. A Gaussian certificate for these
windows must use their actual signed coefficients or prove another valid
comparison. This limitation is explicitly checked and recorded.

The exact verification enumerates path and cycle boundary assignments,
minimizes over all free-edge assignments for paths and cycles through eight
vertices, checks the cubic endpoint formula against the existing path
optimizer, and checks odd/even truncation bounds through order nine. The
general proof above is independent of the tested size cutoff.

## 5. Centering and an explicit bounded-path endpoint charge

Although `H_K` is not generally submodular, the full optimized split cost
`F_G` of a finite path or cycle is submodular and complement-symmetric.
These properties follow from the same partial-minimization argument as in
`PROOF.md`. They give a useful centering lemma.

**Centering lemma.** If `X` is any centered Gaussian vector with unit
coordinate variances and `F` is submodular with `F(b)=F(1-b)`, then

\[
\mathbb E F(\mathbf1[X\le t])\le\mathbb E F(\mathbf1[X\le0])
\qquad\text{for every real }t.
\]

For a positive-definite covariance, differentiation in the common threshold
expresses the derivative as the sum over `i` of

\[
\varphi(t)\,
\mathbb E[\Delta_iF(\mathbf1[X_{-i}\le t])\mid X_i=t].
\]

The conditional Gaussian residual distribution is independent of `t`, and
its effective coordinate thresholds are `(1-rho_ij)*t`. These thresholds
are nonnegative for `t>=0`. The discrete derivative `Delta_iF` is decreasing
in every other bit by submodularity. At `t=0` its conditional expectation
is zero, because complement symmetry gives
`Delta_iF(1-b)=-Delta_iF(b)` and the conditional residual is symmetric.
Thus the derivative is nonpositive for `t>=0`. Complement symmetry treats
negative thresholds. Singular covariances follow by unit-diagonal
positive-definite approximation and bounded convergence; each coordinate
has a continuous marginal. This proves the lemma.

Consequently a centered Gaussian window estimate can bound `F_G` at zero,
after which the centering lemma controls its cost at all thresholds. It is
unnecessary to assert that the non-submodular `H_K` itself is maximized at
zero.

Now consider a cubic path on `k>=K+1` vertices, with odd `K`. Let `b_i` be
the boundary bit supplied by the original perfect matching at vertex `i`.
The two additional endpoint boundary bits can force at most two split
vertices. The first `K-1` path edges contain a matching of size at most
`(K-1)/2`. Apply the odd-window bound to every remaining edge. The result is
the exact deterministic inequality

\[
F_{P_k}\le\frac{K+3}{2}
+\sum_{j=0}^{k-K-1}H_K(b_j,\ldots,b_{j+K}).               \tag{1}
\]

There are `k-K` complete windows. Thus a uniform centered expectation bound
`E H_K<=beta` gives

\[
\mathbb E F_{P_k}\le(k-K)\beta+\frac{K+3}{2}
\le k\beta+\frac{K+3}{2},\qquad\beta\ge0.                \tag{2}
\]

No independence of the endpoint or matching scores is assumed.

## 6. A conditional finite-window reduction to kernel cutwidth

Fix odd `K` and nonnegative `beta`. The following two uniform Gaussian
statements would suffice:

1. Every window of `K+1` consecutive distinct vertices on a 2-factor cycle,
   with its original perfect-matching edge scores, satisfies
   `E H_K<=beta` at threshold zero. Its admissible covariances must include
   matching chords and repeated boundary keys.
2. Every 2-factor cycle of length at most `K` has centered total optimized
   split cost at most `beta` times its number of vertices. For the required
   simple triangle-free kernels, only lengths `4,...,K` occur.

The uniform quantifier is over all required kernels, chosen perfect
matchings, and indicated windows or cycles. These statements may be proved
directly for the finite-radius Gaussian fields used in `PROOF.md`, or for
their admissible cone covariances with a uniform vanishing shell error.
More explicitly, it suffices that the finite-radius fields satisfy the
displayed bounds with `beta+epsilon_R`, where `epsilon_R` tends to zero
uniformly as the field radius `R` tends to infinity. The formulas below
suppress that common vanishing addition.

For the cone formulation, the transfer must use a uniform mean-square
coupling of the finite-radius scores with the cone scores and the resulting
probability of a sign change in one of the finitely many coordinates.
Boundedness of `H_K` then controls its expectation error. Coordinatewise
Gaussian monotonicity cannot be used for `H_K`, which is not submodular.
Threshold centering is applied only to the full optimized gadget cost `F_G`.

The first condition covers longer cycles: each complete cyclic window has
distinct core vertices, so the exact cyclic inequality gives total cost at
most `beta` times cycle length. The second condition handles the shorter
cycles without treating repeated positions as distinct Gaussian variables.

To retain bounded gadgets, choose an integer `L>=K+1`. Keep every cycle of
length less than `L` as one whole-cycle gadget. Partition every longer cycle
into contiguous paths whose lengths lie in `[L,2L-1]`, keeping the cut cycle
edges as base edges. Such a partition exists: if the cycle has length `n`,
take `q=floor(n/L)` parts with sizes differing by at most one. The path
endpoint charge in (2) is at most `(K+3)/2` per part. Since every such part
has at least `L` vertices, the total centered expected cost is at most

\[
\left(\beta+\frac{K+3}{2L}\right)|V(H)|.                 \tag{3}
\]

Every gadget has at most `2L-1` vertices, with `L` fixed before graph order
tends to infinity. The centering lemma and the existing bounded-gadget
insertion theorem therefore apply, including monotone least choices, bounded
tie changes, the vanishing field-shell error, and concentration. For each
fixed `L`,

\[
\operatorname{cw}(H)
\le\left(\beta+\frac{K+3}{2L}+o_L(1)\right)|V(H)|.
\]

For any fixed slack, choose `L` large enough, then a field radius `R` large
enough for the uniform shell error, and finally let graph order tend to
infinity so that the bounded-dependency concentration error vanishes.
This yields the universal asymptotic coefficient
`kappa<=beta`, conditional on the two Gaussian statements above. No
unbounded-cycle concentration claim is required.

## 7. A conditional certificate with a telescoping correction

Fix a real-valued function `Q` on `K` ordered boundary bits and define

\[
W_Q(b_0,\ldots,b_K)
=H_K(b_0,\ldots,b_K)
+Q(b_0,\ldots,b_{K-1})-Q(b_1,\ldots,b_K).
\]

Use this same function, with a consistent orientation, at every window
around a cycle. Each suffix is the next window's prefix. Therefore the
potential terms cancel pointwise, and

\[
\sum_{i\ \mathrm{cyclic}} W_Q(b_i,\ldots,b_{i+K})
=\sum_{i\ \mathrm{cyclic}} H_K(b_i,\ldots,b_{i+K}).
\]

This identity also holds when score coordinates repeat: it uses the actual
ordered bits, without an independence assumption. Thus a uniform centered
bound `E W_Q<=beta` can replace `E H_K<=beta` in the first Gaussian
hypothesis of Section 6.

For a path with `k>=K+1` vertices, there are `k-K` complete windows and

\[
\begin{aligned}
\sum_{j=0}^{k-K-1} H_K(b_j,\ldots,b_{j+K})
={}&\sum_{j=0}^{k-K-1} W_Q(b_j,\ldots,b_{j+K})\\
&-Q(b_0,\ldots,b_{K-1})
 +Q(b_{k-K},\ldots,b_{k-1}).
\end{aligned}
\]

Write `osc(Q)=max Q-min Q`. The last two terms are at most `osc(Q)`, so
the exact path endpoint charge in (1) increases by at most `osc(Q)`.
Consequently the same partition into paths of length at least `L` gives

\[
\operatorname{cw}(H)
\le\left(\beta+
 \frac{(K+3)/2+\operatorname{osc}(Q)}{L}
 +o_L(1)\right)|V(H)|.
\]

With `Q` fixed, taking `L` large gives `kappa<=beta` under the corrected
window hypothesis and the unchanged whole-cycle hypotheses for lengths
at most `K`. The uniform finite-radius or limiting-cone conditions of
Section 6 still apply. Since `Q` is a finite table, its range is bounded,
so the same sign-change coupling controls its shell error. Threshold
centering is still applied to the full submodular cost, not to `W_Q`.

This interface may overcome an obstruction at an isolated relaxed window:
the correction can charge incompatible prefix and suffix distributions
while cancelling on every actual cycle. It cannot lower the sum for a
genuine cyclic configuration, and any proposed `Q` must satisfy one uniform
inequality over the full required domain. No useful potential or corrected
Gaussian bound is claimed here.

## 8. Frontier comparison and a proved window obstruction

For comparison, consider the hypothetical target

\[
K=7,\qquad\beta=\frac{67}{500}=0.134.
\]

It would require an eight-score window bound and short-cycle bounds for
lengths `4,5,6,7`. If proved on all actual score fields, with the uniform
finite-radius control specified above, the circuit conversion would give
`c>=317/67`, exceeding the currently reported `6341/1341`.

**The bare uniform window target is false on the relaxed P8 tree face.**
[WINDOW_OBSTRUCTION.md](WINDOW_OBSTRUCTION.md) gives an exact positive
rational mass vector for which the expectations of forward H7, reversed
H7, and their average each strictly exceed

\[
\frac{17}{125}=0.136>0.1341.
\]

The certificate only needs the pointwise lower bounds
`H7>=D_6(1-D_5)` and `H7_reversed>=D_0(1-D_1)`. Their three-score Gaussian
expectations are enclosed by interval evaluation of the pair-sign arcsine
formula; each is greater than `0.136076509488864`. No high-dimensional
integration is needed. The earlier independent numerical searches and
holdouts are recorded separately in
[WINDOW_ADVERSARIAL.md](WINDOW_ADVERSARIAL.md).

The conditional theorem above remains valid for actual score fields, but
its bare hypothesis cannot be certified over this whole relaxed face.
An abstract distance-type mass vector need not be realized by an actual
graph, and the certificate does not bound the spatial average on a cycle.
A viable stronger result
may need restrictions proved from graph realizability, an inequality for
the average around a cycle, or a useful potential in Section 7. No such
Gaussian bound is established here. The exact
window identities and bounded-gadget reduction specify the Boolean and
graph components to retain when changing that Gaussian proof obligation.

## Replay

```sh
python -B research/gadget_insertion/cycle_windows.py \
  --output research/gadget_insertion/cycle_window_checks.json
```
