# Bounded-gadget insertion: analytic reduction and outstanding certificates

## Status and scope

This note proves an extension of the **insertion mechanism** in `proof.md`,
Sections 2.5 and 3.2–3.4. It specifies finite Gaussian inequalities which would
be sufficient for a stronger universal cubic cutwidth bound. **No new numerical
upper bound for all cubic kernels is proved here.** Tree evaluations, sampled
cone points, and floating-point optimization are exploratory evidence only.

The Gaussian threshold and correlation comparison lemma used below is Lemma G
of `proof.md`. The circuit-to-cutwidth implication is Theorem K2C there. Their
hypotheses and normalization are retained: a bound

    cw(H) <= (kappa + o(1)) |V(H)|

uniformly for the simple, triangle-free, 2-connected cubic kernels implies

    liminf C_B2(f_n)/n >= 1 + 1/(2 kappa).

The upstream notation is `M = 2 kappa`; a per-vertex gadget cost is a candidate
for `kappa`, not for `M`.

For the circuit application there is a particularly concrete sufficient
target: certify the four path inequalities for **P_4, P_5, P_6, and P_7** in
Section 5. The kernel is triangle-free, so choosing L=4 covers it using these
four shapes alone. No short-cycle gadget is needed for this specialization.

## 1. A general gadget and its Boolean cost

Let H be a finite simple cubic graph. Partition its vertices into nonempty
blocks B, each of size at most K. For each block choose a set I_B of *free*
edges with both endpoints in B. Different blocks have disjoint free edges.
All remaining edges form the set E_0 of *base* edges. These choices depend on
H, but must be made independently of the Gaussian noise used subsequently.

For each B, the boundary coordinates may be indexed either by the distinct
base edges incident to B or by their incidences at vertices of B. The incidence
convention repeats a coordinate when a base edge has both endpoints in B.
In an actual assignment the two occurrences have the same bit. We first use
distinct base-edge coordinates.

For a free-edge assignment x and a boundary assignment s, define

    g_B(x,s) = sum_{v in B} 1[the three edge bits incident to v are not all equal],
    F_B(s)   = min_{x in {0,1}^{I_B}} g_B(x,s).

There are no free variables when I_B is empty; the minimum then means the
single available value. All costs are integers between 0 and |B|.

We require the following endpoint condition: each connected component of
(B,I_B) that contains a free edge has a vertex incident to a base edge. This
condition holds in the path/cycle construction in Section 5, because every
vertex has an incident edge of a perfect matching which is kept in E_0.
Without this condition an unanchored component can stay on side zero even
when all boundary bits are one, so its entry time need not be finite.

### Lemma 1: submodularity, symmetry, and monotone least minimizers

The function g_B is submodular in all its free and boundary bits and invariant
under simultaneous complementation of all bits. The partially minimized
function F_B is also submodular and complement-invariant. For each boundary
assignment s, the set of minimizing free assignments has a coordinatewise
least element x_*(s), and

    s <= t  implies  x_*(s) <= x_*(t).

Under the endpoint condition, x_*(0)=0 and x_*(1)=1.

**Proof.** For three Boolean bits,

    1[a,b,c not all equal] = (|a-b| + |a-c| + |b-c|)/2.

Each term |a-b| is a submodular cut function. Their sum proves submodularity
of g_B; complementation preserves every summand. Identifying several
coordinates also preserves submodularity: the coordinate-repetition map
preserves componentwise meet and join. This handles any repeated boundary
occurrences without an independence assumption.

For minimizers x of g_B(.,s) and y of g_B(.,t), submodularity gives

    F_B(s) + F_B(t)
      >= g_B(x AND y, s AND t) + g_B(x OR y, s OR t)
      >= F_B(s AND t) + F_B(s OR t).

This is submodularity of F_B. Complement x as well as s to obtain symmetry.
With s=t, the same inequality shows that both x AND y and x OR y are
minimizers. The finite set of minimizers is therefore closed under meet and
has a least element, its componentwise meet.

If s<=t, take x=x_*(s) and y=x_*(t). The first inequality becomes

    F_B(s)+F_B(t) >= g_B(x AND y,s)+g_B(x OR y,t) >= F_B(s)+F_B(t).

Consequently x AND y is a minimizer at s. Leastness of x gives x<=x AND y,
and hence x<=y.

The all-zero assignment has cost zero and is the least possible assignment,
so x_*(0)=0. At all-one boundary, the all-one free assignment has cost zero.
In any zero-cost assignment, a vertex incident to a base edge of bit one
must have all its incident bits one. This property propagates along the
free edges of its connected component. The endpoint condition therefore
forces every free bit to be one, proving x_*(1)=1. QED.

This proof also supplies an exact min-cut formulation: use a vertex for each
edge-bit variable and a weight-1/2 edge for each pair appearing at an original
cubic vertex, with fixed boundary bits clamped to their prescribed sides.
Parallel pair contributions are added. Exact enumeration or an exact min-cut
routine can verify F_B; no Gaussian computation is involved in this lemma.

## 2. Entry times, ties, bags, and cutwidth

Give each base edge b a finite real score k_b, with ties allowed. Let

    s_B(theta)_b = 1[k_b < theta].

For a free edge e in I_B, set

    tau_e = inf {theta : x_*(s_B(theta))_e = 1};

for a base edge set tau_b=k_b. By Lemma 1, every free coordinate is a monotone
step function which starts at zero and ends at one. It changes only when a
boundary score is crossed. In particular, tau_e is one of the finitely many
base scores incident to B, and the strict convention gives the exact identity

    1[tau_e < theta] = x_*(s_B(theta))_e

for every real theta, including theta equal to a score. Several free or base
edges may enter at the same time; this causes no difficulty.

Order all edges by (tau_e, a fixed distinct edge index). Write first(v) and
last(v) for the first and last positions of edges incident to v, and define

    X_i = {v : first(v) <= i <= last(v)},
    Cost(theta) = #{v : v has an edge with tau < theta and one with tau >= theta}.

Every edge at position i has both endpoints in X_i, and the bags containing
each vertex form an interval. Thus the X_i form a path decomposition. Since
the blocks partition the vertices and the optimized assignment is realized
by the entry times,

    Cost(theta) = sum_B F_B(s_B(theta)).                                      (1)

### A uniform bound on tied-score and grid errors

Let m=max_B |I_B|. For any set A of real numbers, each free edge whose entry
time lies in A can be charged to an incident base edge whose score is that
entry time. A base edge is incident to at most two blocks, and each block
has at most m free edges. Therefore

    #{e : tau_e in A} <= (1+2m) #{b in E_0 : k_b in A}.                       (2)

This remains valid when a base edge has both endpoints in the same block:
it is then incident to just one block. It also remains valid under arbitrary
ties. Define D_0=2(1+2m). For a grid theta_0<...<theta_J, let W_j count base
scores in [theta_j,theta_{j+1}), W_- those below theta_0, and W_+ those at
least theta_J. Every bag satisfies one of

    |X_i| <= D_0 W_-,
    |X_i| <= D_0 W_+,
    |X_i| <= Cost(theta_j) + D_0 W_j for some j.                              (3)

To prove the third alternative, let tau of the edge at position i lie in
[theta_j,theta_{j+1}). If v in X_i has an edge with entry time below theta_j,
it is counted by Cost(theta_j). Otherwise it has an incident edge entering
in that window. Bound the number of such endpoints using (2). The tails
follow by the same argument. No assumption that scores or entry times are
distinct was made.

### Cubic cutwidth conversion

Order the vertices by the position of their median incident edge, breaking
ties by a fixed vertex index. For any nonempty proper prefix, let t be the
median position of its last vertex. A crossing edge at a position below t is
the unique lowest edge of its endpoint outside the prefix; charge it to
that endpoint, which belongs to X_t. A crossing edge above t is the unique
highest edge of its endpoint inside the prefix; charge it there, also in
X_t. These two kinds of recipients lie on opposite sides, and no vertex
receives two charges. There is at most one edge at position t. Consequently

    cw(H) <= max_i |X_i| + 1.                                                (4)

This is the argument of upstream Lemma A.6(e), applied to the new edge order.
It uses cubic degree, not a claimed equality between arbitrary graph layout
parameters.

## 3. Gaussian expectation and concentration: the conditional theorem

Use exactly the upstream radius-R field and normalized edge keys:

    a_v(z) = 2^{-dist(v,z)/2} 1[dist(v,z)<=R],
    a_hat_v = a_v / ||a_v||,
    y_uv = (a_hat_u+a_hat_v)/||a_hat_u+a_hat_v||,
    k_uv = <y_uv,omega>,  omega_z independent N(0,1).

Every base score is standard normal. They need not be independent and may
coincide identically. Let C_B be the covariance of the boundary-score
occurrences of B. By Lemma 1 and upstream Lemma G(a),

    E F_B(s_B(theta)) <= Phi_B(C_B),
    Phi_B(C) := E_{Z~N(0,C)} F_B(1[Z<0]).                                   (5)

The lemma allows singular correlation matrices, so repeated boundary scores
are covered. The objective on the right uses the incidence-indexed version
of F_B with no imposed equalities; the actual singular Gaussian enforces
the required equalities automatically.

**Conditional insertion theorem.** Fix K and a finite collection of connected
gadget shapes, with the endpoint condition above. Suppose every kernel H
has a noise-independent partition into such gadgets. Suppose there are
nonnegative bounds beta_B, depending only on the gadget shape, such that
uniformly over all applicable H and all its gadgets,

    Phi_B(C_B) <= beta_B + epsilon_R |B|,       epsilon_R -> 0,
    sum_B beta_B <= kappa |V(H)|.                                           (6)

Then for every xi>0 there is h_0(xi) such that every such kernel with h>=h_0
has cutwidth at most (kappa+xi)h+1. In particular the remainder is uniform
over the graph class.

**Proof.** First choose R large enough to make epsilon_R small compared with
xi. The contribution F_B lies in [0,K] and depends only on noise within
distance R+K of a selected vertex of B, because the gadget is connected and
each boundary endpoint is within distance K. Distinct blocks have distinct
selected vertices. Contributions whose selected vertices are more than
2(R+K) apart are independent. A ball of fixed radius in a cubic graph has
bounded size, so Var Cost(theta)<=C_{R,K}h uniformly in theta and H. The same
variance bound holds for each base-score window or tail count: assign each
base edge to one endpoint and group its indicator there; each grouped term
lies in [0,3] and has a bounded dependency radius.

By (1), (5), and (6), E Cost(theta)<=(kappa+epsilon_R)h for every theta.
Since |E_0|<=3h/2 and every score is standard normal, a window of width
Delta has mean count at most (3h/2)Delta/sqrt(2pi). Choose Delta small and
T large so that D_0 times every window mean and each tail mean is small
compared with xi h. Use a finite grid covering [-T,T]. Its number of points
depends on xi and K, never on h. Chebyshev and a union bound show that, for
all sufficiently large h, there is one noise realization in which all grid
costs and all window/tail counts differ from their means by less than a
sufficiently small fixed multiple of xi h. Apply (3) and (4).

The quantifier order is: fixed gadget collection and kappa; then xi; then R,
T, Delta, and the finite grid; then h_0; then H, its deterministic partition,
and finally the good noise realization. QED.

Choosing a partition adaptively from the sampled noise is not justified by
this proof. A separate argument would be required for that variant.

## 4. A finite distance-type cone sufficient for the path certificates

For a k-vertex path, k>=2, use core labels 0,...,k-1 and k+2 separate boundary
leaf labels, with owners

    (0,0,1,2,...,k-2,k-1,k-1).

The required edges are the k-1 consecutive path edges and each owner–leaf
edge. This extended labeled graph is a tree with 2k+2 vertices. Labels denote
**occurrences**, not necessarily distinct vertices in H. Core vertices are
distinct. Boundary leaves can coincide with one another or with other core
vertices, subject to the actual simple graph. Thus a non-induced path is
allowed, and a chord not chosen as a free path edge remains a base edge.

For a source z in H, put t_i=dist(pi(i),z)-dist(pi(0),z), where pi maps labels
to actual vertices. These integers satisfy t_0=0 and |t_i-t_j|<=1 on every
required edge. Each core vertex has all three of its actual neighbors
represented in the extended graph. If z is not a core vertex, then

    min_{j in neighbors(i)} t_j = t_i-1              for every core i.       (7)

If z is core s, then all its three neighbors have label t_s+1, every other
core has label strictly greater than t_s and satisfies (7), and every leaf
has label at least t_s. A boundary leaf attached to an owner equal or
adjacent to s along the path cannot itself equal s: that would be a loop or
duplicate path edge. For these owners its label is at least t_s+1.

Let T_k contain every integer pattern meeting the generic or special rules
just described. There are finitely many: enumerate the {-1,0,1} increments
on the extended tree. These necessary conditions may admit unrealizable
patterns. That enlarges the cone and is safe for an **upper** certificate.
The routine `path_cone.enumerate_types` implements this definition; its k=2
cardinality check against 123 is a consistency check, not a replacement for
the source-to-type argument above.

For t in T_k define the positive vector g_t(i)=2^{-t_i/2}, and for nonzero
nonnegative masses m define

    G(m) = sum_t m_t g_t g_t^T,
    G_hat(m) = corr(G(m)),
    C_k(m) = corr(A G_hat(m) A^T),                                         (8)

where A has one row for each boundary occurrence and has entries one at
its core owner and its leaf. Here corr divides an entry by the square root
of its two diagonal entries. Every diagonal is positive; the edge-sum
diagonals before their final normalization are at least two. The matrix
C_k(m) is positive semidefinite with unit diagonal, including singular cases.

### Truncation shell and completeness of the cone relaxation

Choose core 0 as the root. Every extended label lies within distance
D=k of it, and

    A_k := sum_i 2^{dist_extended(0,i)} = 2^{k+2}-1.

For R>k split the actual labeled row Gram matrix into sources at root
distance at most R-k and the remaining sources:

    G_actual = G_in + S.

For a deep source every row is untruncated. Its row contribution is
2^{-dist(root,z)/2}g_{t(z)}, so

    G_in = sum_{t in T_k} m_t g_t g_t^T,
    m_t = sum_{z deep, t(z)=t} 2^{-dist(root,z)}.

The root source is included and has positive mass. Write N_j for the number
of sources at root distance j and s_j=N_j 2^{-j}. In a cubic graph,
s_{j+1}<=s_j for j>=1. A nonzero shell source has root distance between
R-k+1 and R+k. Its squared row norm is at most A_k 2^{-j}, even if some
labels coincide. Consequently

    tr S <= 2k A_k s_{R-k},
    tr G_in >= sum_{j=0}^{R-k} s_j >= (R-k)s_{R-k},

and hence

    tr S <= [2k(2^{k+2}-1)/(R-k)] tr G_in.                                 (9)

For k=2 this is exactly the upstream factor 60/(R-2).
Since |t_i|<=dist_extended(0,i)<=k,

    min_i (G_in)_{ii} >= tr G_in / [2^k(2^{k+2}-1)].                       (10)

Equations (9)–(10) show that the shell is a uniformly vanishing relative
perturbation of every diagonal. The two normalizations in (8) are uniformly
continuous on these matrices after scaling their trace to one. Therefore
the actual boundary covariance differs from C_k(m) by o_R(1), uniformly over
H and the path occurrence. This can also be proved by the entrywise estimates
of upstream Lemma A.3 with the diagonal constant (10).

For fixed dimension, Phi_B is uniformly continuous on the compact set of
all correlation matrices. Indeed, convergent covariance matrices give weakly
convergent Gaussian distributions; the bounded sign-cost function is
continuous off coordinate-zero hyperplanes, each of probability zero since
all variances are one. This is valid even at singular matrices. The shell
therefore changes the expectation by o_R(1).

The following **unproved finite-dimensional inequality** would now suffice
for a per-path bound beta_k:

    Phi_{P_k}(C_k(m)) <= beta_k      for every m_t>=0 with sum_t m_t=1.     (11)

Homogeneity permits the normalization sum m_t=1. It is essential that (11)
cover the entire simplex, including its boundary and singular covariances.
A value at the tree distribution is only one feasible value of its left
side, and thus gives a lower bound on this supremum, not an upper bound.

## 5. Covering every cubic kernel: paths and short cycles

Fix an integer L>=2. A kernel is bridgeless and cubic, so it has a perfect
matching M by the same Petersen theorem used upstream. Choose such a matching
deterministically. The graph H-M is a spanning disjoint union of simple
cycles. In a triangle-free kernel their lengths are at least four.

For a cycle with n>=L vertices, write n=qL+r with q>=1 and 0<=r<L. Cut the
cycle into consecutive paths with vertex counts L+r,L,...,L (q paths total).
Every resulting size lies in [L,2L-1]. Its free edges are its consecutive
path edges. If q=1, omit the closing cycle edge from the free path and keep
that edge in E_0 even though both its endpoints belong to the same block.

For a cycle with n<L, use the whole cycle as one gadget and make its cycle
edges free. Always keep every edge of M in E_0, including a matching chord
whose two endpoints lie in one gadget. Keep the cut cycle edges in E_0 too.
Every vertex then belongs to exactly one connected block of size at most
K=2L-1, and every vertex is incident to a base edge. Thus the endpoint
condition and all the deterministic insertion arguments apply. One can use
m<=2L-2 in (2) for this construction.

The finite family of shapes that must be certified is therefore

    paths P_k,   L <= k <= 2L-1;
    cycles C_c,  4 <= c < L.                                               (12)

### The concrete circuit target: L=4 and only P_4 through P_7

The kernels of upstream Theorem K2C are simple and triangle-free. Every
cycle of H-M therefore has length at least four. With L=4 the cycle family
in (12) is empty: partition each cycle into paths on 4, 5, 6, or 7 vertices.
For a cycle of length 4 through 7, use one path and retain its closing edge
as a base edge with two boundary occurrences. This is already covered by
the covariance and tie arguments above.

It is thus sufficient, for one positive kappa_*, to prove exactly

    max_{k in {4,5,6,7}} (1/k) sup_{m>=0, sum m=1}
        Phi_{P_k}(C_k(m)) <= kappa_*.                                     (12a)

The resulting kernel cutwidth bound is (kappa_*+o(1))h and the resulting
circuit coefficient is 1+1/(2 kappa_*). This route needs neither a theorem
for cubic graphs with bridges nor a separate cycle-gadget certificate.
The cones described above cover all simple cubic neighborhoods, which is
a sufficient outer relaxation of the narrower kernel class. A P_4 bound
alone does not discharge the P_5, P_6, and P_7 obligations.

### Optional short-cycle family for other choices of L

For each short cycle use c core vertices with their cyclic required edges
and one labeled boundary leaf per core vertex. The same generic/special
source rules and the required-edge Lipschitz conditions give a finite cone.
Use the same covariance construction (8). For its shell proof one may take

    D = floor(c/2)+1,
    A = sum_labels 2^{dist_required(root,label)},
    tr S <= [2D A/(R-D)] tr G_in,
    min diag G_in >= tr G_in/(2^D A).

The cycle Boolean cost minimizes over all c cycle-edge bits. It satisfies
Lemma 1. Repeated matching-edge coordinates, such as a chord whose ends
both belong to the cycle, are again represented by singular covariance.

It is sufficient to certify, for every shape in (12),

    sup_m Phi_{P_k}(C_k(m)) <= kappa k,
    sup_m Phi_{C_c}(C_c(m)) <= kappa c.                                    (13)

Then the sum in (6) is at most kappa h, and the conditional insertion theorem
proves the desired universal bound. More flexible shape-dependent budgets
are possible if one separately proves that the chosen cover has total
budget at most kappa h. Short cycles cannot simply be discarded as a
vanishing fraction: a graph may have linearly many vertices in them for
the selected perfect matching.

### A possible rigorous interface for global averaging

There is a more flexible sufficient condition than bounding every path by
kappa times its size. Orient each cycle of the 2-factor consistently. Each
path block then has an entering cut edge ell_B and a leaving cut edge r_B,
both oriented along that cycle. Suppose P is one consistently defined
deterministic real function of an oriented cut edge and graph/field data.
If a proof establishes

    Phi_B(C_B) <= kappa |B| + P(ell_B) - P(r_B) + epsilon_R |B|

for every block, the correction terms cancel exactly when summed around
each cycle. For a cycle covered by a single path, ell_B and r_B are the same
oriented edge and already cancel within that block. Thus the total expected
cost still has the bound needed in Section 3. A finite cone implementation
would have to include all labels needed to evaluate the two corrections
and prove the inequality uniformly, together with its shell-error estimate.
This is a sufficient certificate design, not a claim that a useful P has
been found. In particular, neither averaging covariance matrices nor
averaging type masses may be substituted without proof for averaging the
nonlinear expectations Phi_B.

### Boundary identifications are not extra free variables

When a remaining edge has both endpoints in a block, its two boundary
occurrences must retain one shared base score. The finite cone relaxation
allows more abstract occurrence covariances than are realizable, but every
actual repeated edge has identical normalized key rows. Thus certifying
(11) or (13) on the full relaxed cone is sufficient. Treating the two sides
of the same edge as independent adjustable free bits would describe a
different, unrealizable ordering and is not justified.

An alternative is to make additional internal matching chords free. Sections
1–3 still apply when all relevant internal components remain anchored, but
the gadget cost, free-edge count, source-type collection, and covariance
certificate must then be changed to include every such topology. The path
code in this directory does not implement that alternative.

## 6. What remains before any improved constant can be claimed

1. Choose a complete covering scheme such as (12), or prove another scheme
   covering every kernel with controlled budgets. The concrete choice
   L=4 already gives the complete path-only scheme (12a).
2. Supply exact Boolean costs for all selected shapes. Current path code
   covers paths. Short-cycle implementation/certification is needed only
   when a different scheme actually includes those shapes.
3. Prove rigorous upper bounds over **all** masses in the corresponding
   finite cones, with outward-rounded arithmetic or another valid analytic
   certificate. Floating-point optimizers and sampled points do not do this.
4. If reducing the cones or using a purported extremizer, prove that the
   reduction covers coincidences, non-induced paths, zero masses, singular
   covariances, and all relevant graph topologies.
5. Retain exact or rigorously enclosed arithmetic when translating a proved
   kappa into the circuit coefficient 1+1/(2 kappa).

The combinatorial insertion mechanism and the sufficient certificate
interface are established above. The new universal numerical certificate
needed for a stronger circuit lower bound remains open in this work.

### Exploratory warning: do not assume the path tree point is extremal

The accompanying exploration has found a nonuniform perturbation of the
P_4 tree-type masses whose numerical expected cost exceeds the value at the
uniformly rooted tree distribution. The final 128-point quadrature reports
per-vertex values about 0.1361262327372273 at the tree distribution and
0.13613542578979634 at the numerical stationary candidate; see
`cone_search_P4.json`. The 64-point values agree to about 10^{-15}, but
agreement between quadrature orders is not an error certificate. These are
floating-point estimates, not certified inequalities. Their role
is to warn against assuming that the upstream tree-extremality argument
extends automatically to larger paths. A small numerical stationarity
residual does not establish a local or global maximum. Even the tree-supported
face needs analysis. A successful certificate might bound a nonuniform extremizer or
use additional realizability or global averaging constraints; those options
require their own proofs.
