# An explicit, very small uniform P_4 improvement

**Historical first quantitative proof.** This note preserves the original
10^(-54) argument. The stronger current result and its revised constants
are in `STABILITY_REFINEMENT.md`; the stronger graph conversion is in
`GLOBAL_REDUCTION.md`. See `improved_bound.json` for their exact numerical
consequence.

## Statement and certificate inputs

Let M=M_tree be the exact upstream pair constant, and let M_4 be the maximum
total P_4 cost over the entire 4,655-type simplex. The arguments below give

    M_4 <= 2M - eta,       eta = 10^(-54).                                 (1)

This bound is deliberately much weaker than the numerical conjecture
M_4/4 approximately 0.13613543. It gives an explicit improvement of the
upstream insertion constant, not a competitive new decimal record.

The proof uses the following named certificate inputs, in addition to the
upstream Gaussian comparison and pair theorem:

1. Upstream R1 witnesses and reductions: every relevant pair covariance C
   either has cost at most Mt_lo, or dominates a zero-determinant witness
   C', which in turn dominates its tight-cut chart matrix C_0. The chart
   has a certified bound

       Phi(C_0) <= M - (3349/10000)||theta-theta_T||^2

   on its local box of radius 10^(-4). Outside the local leaves, the R1
   cover bounds Phi(C_0) by Mt_lo. For the few round-to-nearest operations
   discussed in upstream Section 6.4, we allow the larger upper bound
   Mt_lo+10^(-13), exceeding the stated total error of less than 10^(-14).
2. `research/mtree_refinement/certificate.json` proves

       M >= 0.276583606882373471022763934877,
       Mt_lo = 0.27658360687864314.

   In particular M-(Mt_lo+10^(-13))>3*10^(-12).
3. `certification/p4_tree.json` certifies that every correlation matrix
   within entrywise off-diagonal distance 2^(-20) of the exact natural
   P_4 tree covariance has total cost at most

       4*(6807/50000) = 0.54456.                                          (2)

   This is a local correlation-box theorem, not merely a point quadrature.
   Section 5 of `certification/README.md` also proves that the full-simplex
   mass ball ||m-m_tree||_1<=2^(-30) maps into this correlation box.
4. The exact six-column tree diagonal matrix, on rows (0,1,2,3,4,8), has
   inverse infinity norm 3. Its entries and the exact finite checks are
   specified in Section 5. No floating-point inverse is assumed.

The rest consists of explicit inequalities. The type projections, equality
structure, matching packing, and notation are given in `GLOBAL_STRATEGY.md`.

## 1. A uniform negative derivative near the old tree covariance

Let p=(1+2sqrt(2))/4, s=(2+sqrt(2))/4, rho=2sqrt(2)/3, and let C_T be the
old 4-by-4 tree correlation matrix. On the unit-diagonal box

    max_{i<j}|C_ij-(C_T)_ij| <= 1/1000,                                   (3)

every matrix is positive definite. Indeed lambda_min(C_T)=1-p>0.04, using
sqrt(2)<1.42. The off-diagonal perturbation has operator norm at most 0.003,
so lambda_min(C)>0.03.

Every two-variable conditional covariance S, obtained as a Schur
complement, satisfies S>=0.03 I. This follows by minimizing the quadratic
form of C over the conditioned coordinates. Its diagonal entries are at
most one. After diagonal normalization, its correlation matrix is still
at least 0.03 I. Hence every conditional correlation r_cond is at most 0.97.

The old cross-entry derivative formula is

    -dPhi/dC_ij
      = arccos(r_cond)/(2pi^2 sqrt(1-C_ij^2)).

Since arccos(r_cond)>=sqrt(2(1-r_cond))>=sqrt(0.06)>0.2, pi^2<10, and the
square root in the denominator is at most one, this derivative magnitude
exceeds 0.01. For a pair entry the upstream formula gives magnitude at
least 1/(2pi)>0.01. Thus throughout (3), for all six coordinates,

    dPhi/dC_ij <= -1/100.                                                  (4)

This step is analytic; no numerical derivative enclosure is needed.

## 2. Near equality of a pair cost forces its covariance near C_T

Suppose an old-cone covariance C has Phi(C)>M-eta. The gap in input 2
excludes the outside-thin-set case and every nonlocal R1 leaf. Apply the
upstream witness and tight-cut reduction, using the same coordinate
permutation throughout, to obtain

    C_0 <= C' <= C,
    M >= Phi(C_0) >= Phi(C') >= Phi(C) > M-eta.

The local quadratic bound implies

    ||theta-theta_T|| < 2sqrt(eta),

because 3349/10000>1/4. Put v=sqrt(eta). The pair entries x,y of C_0
therefore differ from p by less than 2v and its two imbalance parameters
have magnitude less than 2v. In the upstream chart,

    S_cut^2 = K_0-a_ab(x+y),       a_ab<5,

and both S_cut and its tree value exceed 3/2. Thus
|S_cut-S_tree|<7v. Each factor S_cut/2 plus or minus an imbalance differs
from sqrt(s) by less than (11/2)v. For v<=10^(-4), each such factor is
less than 1.1. Since sqrt(s)<1, each cross product differs from s by less
than (11/2)(1.1+1)v<12v.
Consequently

    ||C_0-C_T||_max <= 12sqrt(eta).                                       (5)

Put D=C-C_0, which is entrywise nonnegative off the diagonal. The straight
segment C_0+tD consists of correlation matrices. If max D_ij were at least
1/2000, its initial segment ending at maximum entry increment 1/2000 would
lie in box (3), by (5) and 12sqrt(eta)<1/2000, and would lose at least

    (1/100)(1/2000) = 1/200000

in Phi by (4). Monotonicity on the rest of the segment would contradict
Phi(C_0)-Phi(C)<eta. Therefore max D_ij<1/2000 and the whole segment lies
in (3). Integrating (4) once more gives max D_ij<=100eta. Combining with
(5), and using eta<=10^(-4), proves the convenient bound

    ||C-C_T||_max <= r,       r := 32sqrt(eta).                             (6)

The argument covers a singular original C: before any possible exit the
segment lies in the positive definite box, and Gaussian monotonicity covers
the remaining part. The actual conclusion places C itself in that box.

## 3. Two nearby pair covariances control all row norms

Now take a P_4 cone mass vector m with both old pair costs greater than
M-eta. By (6), both pair covariances are within r of C_T. Let q_i=G_ii,
n_i=sqrt(q_i), and normalize m so that q_0=sum m=1. Put u=sqrt(r).
For the chosen eta, r<=10^(-20) and u<=10^(-10).

For either old pair, let Y be the sum of its four boundary unit keys and
Z its internal unit key. The sibling correlations are within r of p.
The star lower bounds give

    <Z,Y> >= 4p-2r.

The covariance bounds give ||Y||^2<=(4p)^2+12r, hence ||Y||<=4p+2r since
p>3/4. Subtracting the lower bound for the opposite star shows that the
mean of the three incident key correlations at either core vertex is
between p and p+4r/3. Thus this is true at all four P_4 core vertices.

In the star function mu(e_1,e_2,e_3), every partial derivative is at least
1/5: writing c_i=sqrt((1+e_i)/2) and c_0=sqrt((1+rho)/2), the derivative
is (sum c_j)/(12c_i)>=c_0/4>1/5. The star inequality therefore implies

    sum_{edges at a core}(e-rho) <= 20r/3 < 7r.                            (7)

In particular every required edge has e-rho<=7r. For such an edge, let
t=n_j/n_i. The type increments imply 1/sqrt(2)<=t<=sqrt(2). Lemma EC gives

    (t-1)^2/t <= 2(e-rho)/rho < 15r,

using rho>14/15. Since t<3/2, it follows that |t-1|<5u. Every label is at
required-tree distance at most four from root 0. Multiplying at most four
ratios and using u<=10^(-10) gives

    |n_i-1| <= 30u,       |q_i-1| <= 61u.                                 (8)

For example (1+5u)^4<=1+30u and (1-5u)^4>=1-30u for u<=1/1000; the
second bound follows from |n_i^2-1|<=60u+900u^2<=61u on this range.

## 4. Star residuals bound the entire mass outside the six tree types

At a core vertex, let U be its normalized row, U_i its normalized neighbor
rows, and Y_i its three unit edge keys. Put c_i=sqrt((1+e_i)/2),
c_0=sqrt((1+rho)/2), and T=sum Y_i-(sum c_i)U. The previous mean-key bound
and the exact identity 3+6p=9c_0^2 give

    T perpendicular to U,       ||T||^2 <= 8r.

Moreover (7) and c_0>0.98 imply sum(c_i-c_0)<2r. Since

    sum U_i-2sqrt(2)U
      = 2c_0(sum Y_i-3c_0 U)+2 sum(c_i-c_0)Y_i,

its norm is at most 6u+8r<8u. Return to the unnormalized rows a_i=n_i U_i.
By (8), n_i<=2, and adjacent norm differences are at most 10u. Hence at
each core vertex v,

    ||R_v|| <= 200u,
    R_v := sum_{w neighbor of v} a_w-2sqrt(2)a_v.                          (9)

The constant 200 is conservative; the preceding estimates already give 46.

For a supported type t, each neighbor/core ratio belongs to
{1/sqrt(2),1,sqrt(2)}. A sum of three such ratios equals 2sqrt(2) precisely
when one is sqrt(2) and two are 1/sqrt(2). Every other possible sum differs
from 2sqrt(2) by more than 1/10. This is a finite exact check of the ten
count triples; its smallest nonzero absolute difference is
3/sqrt(2)-2>1/10.

Also every core entry g_t(v) is at least 1/4. Therefore every type outside
the six tree types has, at some core vertex, a source-coordinate residual
of magnitude at least 1/40. If q_off denotes the total mass outside those
types, orthogonality of the source coordinates gives

    q_off/1600 <= sum_v ||R_v||^2 <= 4*(200u)^2,
    q_off <= 3*10^8 r.                                                     (10)

There is no cancellation between source types in this squared-norm identity.
The characterization of the six zero-residual types is the oriented-tree
argument in `GLOBAL_STRATEGY.md`, Section 4.

## 5. Exact diagonal inversion controls the full mass vector

Let w be the six tree-type masses, in sorted type order, and w_T their
natural values (1/24,1/24,1/12,1/6,1/3,1/3). Let A be their diagonal
matrix on rows (0,1,2,3,4,8). Its exact inverse is

    [   0,    0,     0, -1/24,    0,  1/12 ]
    [   0,    0, -1/12,  5/24,    0, -1/12 ]
    [   0, -1/6,  5/12,  -1/6,    0,     0 ]
    [ -1/3,  5/6,  -1/3,     0,    0,     0 ]
    [ -1/3,    0,     0,     0,  2/3,     0 ]
    [  5/3, -2/3,     0,     0, -2/3,     0 ].

Thus ||A^(-1)||_infinity=3. Every squared generator entry on a selected
diagonal is at most 16. Equations (8) and (10) consequently imply

    ||w-w_T||_infinity <= 3*(61u+16q_off),
    d := ||m-m_tree||_1 <= 1098u+289q_off <= 2000u,                        (11)

where the last inequality uses r<=10^(-20). These include the mass outside
the tree face and do not silently renormalize the six retained masses.

## 6. Entering the certified natural neighborhood

For eta=10^(-54), equation (6) gives r=32*10^(-27). Thus sqrt(r)<2*10^(-13),
and (11) proves

    ||m-m_tree||_1 <= 2000sqrt(r)
                    < 4*10^(-10) < 2^(-30).                             (12)

The mass-to-correlation estimate in input 3 applies on this full-simplex
ball, including nonzero masses on all 4,655 types. It places C_4(m) within
entrywise off-diagonal distance 2^(-20) of the exact natural covariance.
No restriction to the six-dimensional tree face is used here.

## 7. The global dichotomy and the circuit consequence

For every normalized P_4 mass vector, either one old pair cost is at most
M-eta, or both are greater than M-eta.

In the first case the exact nonnegative-gain decomposition gives
Phi_4<=2M-eta. In the second case, Sections 2–6 place its covariance in the
certified natural box. Input 3 gives Phi_4<=0.54456. The certified lower
bound on M gives

    2M-0.54456 > 43/5000 > eta,

so this case also gives Phi_4<=2M-eta. This proves (1) on the entire closed
simplex, including singular covariance cases.

The P_2/P_4 matching packing in `GLOBAL_STRATEGY.md`, Section 7, therefore
gives the explicit universal asymptotic kernel bound

    cw(H) <= (M/2-eta/7+o(1)) |V(H)|,

and the upstream reduction gives

    liminf C_B2(f_n)/n >= 1+1/(M-2eta/7).                                  (13)

The gain is extremely small. It is not the unproved bound M_4/4<=0.13615,
and it is not a claim to exceed the displayed provisional CircuitBreaker
record. Its role is to turn the qualitative global argument into an
explicit uniform inequality while preserving every graph and asymptotic
quantifier.
