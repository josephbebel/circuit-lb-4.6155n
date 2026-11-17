# A strict universal improvement from four-vertex insertion

## Result and its numerical limitation

This note proves a qualitative global statement about the finite P_4 cone:

    M_4 := max_m Phi_{P_4}(C_4(m)) < 2 M_tree.                              (1)

The maximum is over all normalized nonnegative masses on the 4,655 distance
types defined in `path_cone.py` and `PROOF.md`; it is not restricted to actual
graph neighborhoods. The proof uses the upstream certified pair theorem
`M_123=M_tree`, its recorded equality case, and elementary exact arguments
below. It does not rely on numerical optimization or assume that the
nonuniform stationary candidate maximizes the P_4 cone.

Together with a packing argument, (1) proves that bounded-gadget insertion
gives a universal kernel cutwidth constant strictly smaller than M_tree/2.
The compactness argument in this note does not quantify the gap. The
companion `QUANTITATIVE_GLOBAL.md` supplies the explicit lower bound
2M_tree-M_4>=10^(-54), using the upstream stability certificate and the
separately certified natural P_4 neighborhood. This is a very small uniform
improvement, not the stronger bound M_4/4<=0.13614 or 0.13615 and not a new
leaderboard record. The numerical record-target certificate remains open.

The only large-certificate input is the existing upstream pair theorem.
`verify_rigidity.py` supplies small exact integer/rational checks of the new
type-projection and equality-support algebra; its output is
`rigidity_verification.json`. The argument below explains their meaning.

## 1. Exact pair decomposition and a nonnegative gain

Number the P_4 vertices 0,1,2,3. Its boundary edge scores, in order, are
(a,b,c,d,e,f), with owners (0,0,1,2,3,3). Write x,y,z for its three internal
edge scores or side variables, as appropriate. For Boolean boundary bits s,
define

    A_b(s) = F_2(s_a,s_b,b,s_c) + F_2(b,s_d,s_e,s_f),  b in {0,1}.

Here F_2 is exactly the upstream matched-edge cost. Optimizing the outer
internal bits x and z first gives the pointwise identity

    F_4(s) = min(A_0(s), A_1(s)).                                          (2)

The original normalized Gaussian key of the middle edge is available even
though the new ordering will optimize its side. Let b_y=1[y<0]. The gain is

    Delta(s,b_y) = A_{b_y}(s) - min(A_0(s),A_1(s)) >= 0.                    (3)

Taking expectations in the joint seven-score Gaussian distribution gives

    Phi_4(m) = Phi_L(m) + Phi_R(m) - E_m Delta.                            (4)

The two pair covariances in (4) are legitimate old-cone covariances. To see
this directly, label P_4's core vertices 0,...,3 and boundary leaves 4,...,9.
The left pair uses labels (0,1,4,5,2,6); the right pair uses
(2,3,1,7,8,9). Restrict a distance type to those labels and subtract the new
root coordinate. Every required local edge still has increment -1,0,or 1;
every non-source core vertex has a descending neighbor; and a source in
the pair has exactly the special-source constraints of the old 123-type
cone. A source outside the pair is a generic source for that pair.

For a restriction rooted at label u, a full-cone contribution m_t g_t g_t^T
becomes the projected contribution of mass m_t 2^{-t_u}. Equal projected
types are combined and masses are then normalized. Scaling does not affect
the key covariance. This explains the reweighting needed for the right pair.
The exact script checks all 4,655 projections into the 123-type set on both
sides; each projection in fact reaches all 123 types.

Thus the upstream theorem implies Phi_L<=M_tree and Phi_R<=M_tree. In
particular Phi_4<=2M_tree globally, and the exact identity

    2M_tree-Phi_4(m)
      = (M_tree-Phi_L(m)) + (M_tree-Phi_R(m)) + E_m Delta                  (5)

expresses the gap as three nonnegative terms. We prove that they cannot
all vanish at any cone point.

## 2. Equality in the upstream pair theorem forces its tree covariance

Use the upstream notation

    rho = 2 sqrt(2)/3,
    p = (1+2 sqrt(2))/4,
    s = (2+sqrt(2))/4,
    C_T = [[1,p,s,s],[p,1,s,s],[s,s,1,p],[s,s,p,1]].

The implication needed here is

    Phi_F(C(m))=M_tree  implies  C(m)=C_T.                                (6)

It follows from the written upstream proof, not from an assumed uniqueness
of its mass vector. Outside its thin set N, Proposition NEED bounds the
cost by Mt_lo<M_tree. Inside N, Theorems K-N and PSD-N produce a correlation
matrix C'<=C entrywise with zero cross determinant and the required tree
supporting cut. Theorem R1 states that Phi_F(C')<=M_tree, with equality
only at C'=C_T. Equality in Phi_F(C)<=Phi_F(C')<=M_tree therefore forces
C'=C_T.

If C>=C_T differs in any off-diagonal entry, interpolate from C_T to C.
C_T is positive definite. At C_T every off-diagonal partial derivative of
Phi_F is strictly negative: this follows from upstream equation (2.4),
since its conditional Gaussian correlations lie strictly between -1 and 1.
The initial directional derivative along C-C_T is consequently negative.
Lemma G(b) keeps the function nonincreasing along the rest of the segment.
Hence Phi_F(C)<Phi_F(C_T), a contradiction. This proves (6), including a
possibly singular endpoint C.

## 3. Tree covariance forces exact row and distance-type equalities

Consider one old-cone pair uv with four unit boundary key vectors A,B,C,D
having covariance C_T. Let Z be the unit key vector of uv and put
Y=A+B+C+D. The key correlations within the two boundary pairs are p.
Upstream Lemma STAR gives

    <Z,A>+<Z,B> >= 2p,
    <Z,C>+<Z,D> >= 2p.

Also, exactly,

    ||Y||^2 = 4+4p+8s = (4p)^2.

Cauchy-Schwarz therefore forces <Z,Y>=4p, Z=Y/(4p), and equality in both
star inequalities. The mean of the three incident key correlations at each
core vertex is p.

The chain in Lemma STAR(b) is

    mean key correlation >= mu(e_1,e_2,e_3) >= mu(rho,rho,rho)=p,

where every edge cosine e_i>=rho and mu is strictly increasing in each
argument on this range. Equality of the means therefore forces all five
underlying edge cosines of the pair neighborhood to equal rho.

Lemma EC now forces the two Gram diagonals along each required edge to be
equal: its lower bound is

    e_uv >= (rho/2)(sqrt(G_vv/G_uu)+sqrt(G_uu/G_vv)) >= rho,

with equality in the second inequality only at ratio one. All six local
Gram diagonals are thus equal. Moreover the proof of EC is a sum of
nonnegative per-type terms. Equality forces every supported type to have
increment +1 or -1 on each of the five required edges; a zero increment
would give strictly positive slack.

There is an additional equality in the star projection. Let U_w be the
normalized row at a core vertex and U_1,U_2,U_3 its three normalized neighbor
rows. When all incident edge cosines are rho, equality in STAR gives

    U_1+U_2+U_3 = 3rho U_w = 2sqrt(2) U_w.                                (7)

For completeness, put c_0=sqrt((1+rho)/2). Every incident key has projection
c_0 onto U_w. Equality in the norm-versus-projection step means that their
sum is 3c_0 U_w. Since each key is (U_w+U_i)/(2c_0), rearranging gives (7).

The equal diagonals permit (7) to be read coordinatewise in the source-type
representation. For every type with positive mass, the three neighbor/core
ratios are each either sqrt(2) or 1/sqrt(2), and their sum is 2sqrt(2).
Exactly one ratio is sqrt(2). Equivalently:

    every supported type has exactly one descending neighbor at each core. (8)

All these conclusions concern the outer cone itself; no realizability of
its masses by a finite graph has been assumed.

## 4. Two optimal pair marginals force the natural P_4 tree masses

Suppose Phi_L=Phi_R=M_tree for a P_4 cone point. By Sections 2–3, every
required edge of its extended path has increment +1 or -1 on the support,
and every one of its four core vertices has exactly one descending
neighbor. The two local six-label neighborhoods overlap at core vertices
1 and 2, so their equal-diagonal conditions make all ten Gram diagonals
equal. Normalize their common value to one.

Orient every required edge of a supported type toward the lower label.
The extended path is a tree with four core vertices, six leaves, and nine
edges. Its core outdegrees sum to four. The leaf outdegrees must therefore
sum to five. Exactly one leaf is a sink, and all descending routes lead to
it. Because every step changes the label by one, the type is precisely the
relative distance-to-that-leaf pattern. There are exactly six such types,
the natural tree face.

Let leaf j have owner o_j in (0,0,1,2,3,3), and let its exit-type mass be m_j.
In its own exit type the squared leaf row entry is twice the squared owner
entry. In every other exit type it is half the owner entry. The own-type
squared owner entry is 2^{o_j}, since core 0 is the root. Thus

    G_leaf(j),leaf(j) = (1/2)G_owner(j),owner(j) + (3/2)m_j 2^{o_j}.

Both diagonals are one, forcing

    m_j = 1/(3*2^{o_j}).                                                    (9)

These are the natural tree masses. In sorted type order they are
(1/24,1/24,1/12,1/6,1/3,1/3). The exact diagonal-constraint matrix has rank
six, independently checked by rational elimination in `verify_rigidity.py`.
In particular, simultaneous equality of the two old pair bounds has one
normalized P_4 mass vector: m_tree.

## 5. The natural P_4 point has a strictly positive insertion gain

This step needs no numerical integration.

First the six boundary scores have a nonsingular joint Gaussian law. Put
r=1/sqrt(2), N=sqrt(2+2rho), and let xi_1,...,xi_6 be independent standard
normal variables, one for each exit branch. At the natural masses, the core
field at vertex i is

    V_i = (1/sqrt(3)) sum_j r^{|i-o_j|} xi_j.

The boundary leaf field at j is r V_{o_j}+(r/sqrt(3))xi_j. Consequently the
boundary-score vector has the form

    Z = [(1+r)K+rI] xi / (N sqrt(3)),
    K_ij = r^{|o_i-o_j|}.                                                   (10)

K is positive semidefinite: r^{|i-j|} is the covariance of the elementary
stationary recursion U_{i+1}=rU_i+sqrt(1-r^2)epsilon_{i+1}, and repeating
owner indices preserves positive semidefiniteness. Since r>0, the matrix
(1+r)K+rI is positive definite. Thus Z has positive definite covariance
and every nonempty open box of its six scores has positive probability.

Write alpha=4p=1+2sqrt(2). The tree star equalities imply, for the three
internal edge scores,

    alpha x = a+b+c+y,
    alpha y = x+c+d+z,
    alpha z = y+d+e+f.

Eliminating x and z gives

    y = [a+b+e+f+(alpha+1)(c+d)]/(alpha^2-2).                              (11)

Consider the open event

    a,b in (10,11),       c,d,e,f in (-11/10,-1).                           (12)

The numerator of (11) is greater than

    20-11/5-(11/5)(alpha+1) = (67-22sqrt(2))/5 > 0.

The denominator is positive, so the original middle-edge side is zero.
The boundary bits are (0,0,1,1,1,1). Directly in (2), A_0=2 and A_1=1:
the left pair costs one on either middle side, whereas the right pair costs
one on side zero and zero on side one. Thus Delta=1 on event (12).

By (10), this event has strictly positive probability. Hence E Delta>0 at
m_tree and, by (4), Phi_4(m_tree)<2M_tree.

## 6. Compactness makes the improvement uniform over the entire cone

The normalized mass simplex is compact. Every generator entry is strictly
positive, so all Gram diagonals and edge normalization denominators are
bounded away from zero on this simplex. The map m->C_4(m) is continuous.
The Gaussian sign expectation Phi_4 is continuous even at singular
correlation matrices: every coordinate has variance one and therefore the
coordinate-zero discontinuity hyperplanes have probability zero.
Thus Phi_4 attains its maximum M_4.

If M_4=2M_tree, equality in (5) would make both pair costs equal M_tree and
the gain zero at a maximizing point. Section 4 forces that point to be
m_tree; Section 5 gives strictly positive gain there. This contradiction
proves (1). Define the exact positive constant

    delta = 2M_tree-M_4 > 0.                                               (13)

A useful quantitative certificate interface follows from the same argument.
Let S(m)=(M_tree-Phi_L)+(M_tree-Phi_R) and G(m)=E Delta. The unique zero of
S is m_tree, and G(m_tree)>0. It is sufficient to give an explicit
neighborhood U and positive numbers g,eta such that

    G(m)>=g on U,       S(m)>=eta outside U.                               (14)

Then delta>=min(g,eta). This separates a local gain certificate from an
old-pair stability certificate. Neither positive number is supplied by the
compactness argument alone.

The companion `QUANTITATIVE_GLOBAL.md` uses the same decomposition with a
local bound on the total deficit S+G, rather than on G separately. Its
explicit stability argument gives delta>=10^(-54). That quantitative result
does not require any claim about the numerical maximizer of the P_4 cone.

## 7. A universal graph packing needing only P_2 and P_4

The qualitative theorem does not require P_5,P_6,P_7 certificates.

Let H be a simple bridgeless cubic graph on h vertices and choose a perfect
matching M. Contract every edge of M. The resulting multigraph Q has
q=h/2 vertices and degree four, counting multiplicities. It has no loops:
an additional edge within a matched pair would duplicate the original
matching edge in the simple graph H. Parallel edges between distinct
contracted vertices are allowed.

Choose any maximal matching J in Q and write j=|J|. Its unmatched vertices
U form an independent set. Thus all 4|U| edge incidences from U lead to the
2j matched vertices. Each matched pair of Q-vertices has total degree eight
and at least two incidences used by its chosen matching edge; it can send
at most six incidences to U. Parallel edges only consume more of this
capacity. Therefore

    4(q-2j) <= 6j,       j >= 2q/7 = h/7.                                 (15)

Every edge of J corresponds to an edge of H joining two disjoint edges of
M. Their union is a path on four distinct vertices. These paths are vertex
disjoint. Optimize their three path edges jointly, retaining any additional
chords as base edges. Every unused edge of M is a P_2 gadget. This gives
j P_4 gadgets and h/2-2j P_2 gadgets covering all vertices.

The bounded-gadget theorem in `PROOF.md` applies. The new P_4 cone bound is
M_4; the old pair bound is M_tree. Before the vanishing field-shell and
concentration errors, the total expected cost is at most

    j M_4 + (h/2-2j)M_tree
      = h M_tree/2 - j delta
      <= h (M_tree/2-delta/7).                                            (16)

Consequently every kernel has

    cw(H) <= (kappa_new+o(1))h,
    kappa_new = M_tree/2-delta/7 < M_tree/2.                               (17)

The remainder is uniform over the kernels by the same fixed-gadget
concentration and shell argument as before. This constant is positive:
0<=M_4 implies delta<=2M_tree, hence kappa_new>=3M_tree/14>0.
The upstream circuit reduction yields

    liminf C_B2(f_n)/n >= 1+1/(M_tree-2delta/7) > 1+1/M_tree.               (18)

Equation (18) is a qualitative strict improvement of the upstream method.
It is not a numerical record claim. To approach the current record target,
the path-only cover by P_4 through P_7 and their stronger quantitative
certificates remain an available separate route.

## 8. Why two easier global shortcuts do not suffice

### Optimal old pair marginals alone are insufficient without cone rigidity

For an arbitrary Gaussian PSD gluing of two copies of C_T sharing the
middle score y, the gain can be zero. Let the left remaining triple be
L=r_0 y+epsilon, where r_0=(s,s,p), and use the reflected right triple
R=r_0 y-epsilon. Both pair marginals are C_T. Conditional on y>0, every
entry of R has positive conditional mean and its negative-side bit is at
most the complement of the corresponding left bit. The discrete derivative
of F_2 in its middle side is nonincreasing in its other bits and changes
sign under complementation. Thus the two derivatives sum to a nonnegative
number, and side zero is already optimal. The y<0 case is symmetric.
The gain is identically zero.

This example is a legitimate Gaussian PSD construction. Sections 3–4
exclude it from simultaneous equality in the P_4 distance-type cone.
Any global proof which discards that geometry must replace it with another
valid constraint.

### A single supporting covariance hyperplane is not enough numerically

`saddle_diagnostic.py` computes the 15-coordinate covariance Hessian of
Phi_4/4 at the nonuniform numerical stationary candidate. At difference
steps 10^{-5} and 5*10^{-6}, it finds four positive eigenvalues, approximately

    0.086315642, 0.097548571, 0.256619713, 0.256619713.

Restricting to the gradient-orthogonal tangent plane still leaves four
positive numerical eigenvalues, approximately

    0.086315642, 0.097042114, 0.256619713, 0.256619713.

These are reproducible floating-point curvature diagnostics, not interval
certificates; see `saddle_diagnostic_P4.json`. They explain why a direct
transplant of the upstream one-saddle argument is inadequate. The larger
cone-aware proof remains necessary. Likewise, multistart optimization which
returns the same candidate is evidence for investigation, not a cover of
the full simplex.
