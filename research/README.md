# Cubic cutwidth research: stronger proof, a method obstruction, and longer paths

This work starts from upstream commit
[`a56500014ba159397903dea81d41ac836f6880d5`](https://github.com/rrwilliams0/circuit-lb-4.6155n/tree/a56500014ba159397903dea81d41ac836f6880d5).
All additions are in `research/`; the upstream proof and certificates retain
their original constants.

## Current results

The latest proof establishes

\[
M_4\le2M_{\rm tree}-\eta,\qquad
\eta=\frac{17}{5\cdot10^{19}}=3.4\cdot10^{-19},
\]

and a new graph reduction gives

\[
\kappa\le\frac{M_4}{4}
\le\frac{M_{\rm tree}}2-8.5\cdot10^{-20},
\qquad
c\ge1+\frac1{M_{\rm tree}-\eta/2}.
\]

Here `cw(H) <= (kappa+o(1)) |V(H)|` is the uniform kernel cutwidth bound,
and the upstream circuit reduction is `c=1+1/(2*kappa)`. Smaller `kappa`
means a stronger circuit lower bound. Exact outward decimal consequences are

```
kappa <= 0.138291803441186735426381967438936496879481985
c     >= 4.615543275582792607402913476760178553046176686.
```

**This is a strict global improvement of the original method, with a very
small absolute gain. It does not establish a new leaderboard record.**
The refreshed reported frontier is substantially stronger, and we now have
a rigorous obstruction showing that the present uniform P4 relaxation
cannot reach it. Longer paths have promising natural-tree values, but an
adversarial search led to a second rigorous obstruction: the bare H7 window
bound also fails on its relaxed tree face. The new conditional reduction
allows corrections that cancel across neighboring windows; a successful
corrected Gaussian certificate remains open.

| Result | Status | Main reference |
|---|---|---|
| Explicit P4 deficit `17/(5*10^19)` | Analytic proof with exact finite checks, using named upstream certificates | [STABILITY_REFINEMENT.md](gadget_insertion/STABILITY_REFINEMENT.md) |
| Transfer `kappa <= M_4/4` using P4 alone | Graph proof; handles all 2-edge-cut chains with a uniform additive constant | [GLOBAL_REDUCTION.md](gadget_insertion/GLOBAL_REDUCTION.md) |
| Reduction to cyclically 4-edge-connected cubic graphs | Further graph-class reduction; no numerical coefficient claimed | [GLOBAL_REDUCTION_CYCLIC3.md](gadget_insertion/GLOBAL_REDUCTION_CYCLIC3.md) |
| Natural P4 mass neighborhood of L1 radius `0.0036` | Interval derivative certificate and a wider covariance comparison certificate | [MASS_NEIGHBORHOOD.md](gadget_insertion/certification/MASS_NEIGHBORHOOD.md) |
| Ceiling below the refreshed frontier for the uniform P4 route | Exact rational consequence of an already certified feasible point | [FRONTIER_OBSTRUCTION.md](gadget_insertion/FRONTIER_OBSTRUCTION.md) |
| Joint path/cycle cost as a matching problem; finite-window bounds | Exact Boolean identities with exhaustive finite checks | [CYCLE_WINDOWS.md](gadget_insertion/CYCLE_WINDOWS.md) |
| Obstruction for forward, reverse, and averaged H7 windows | Exact rational covariance witness and interval arcsine bounds | [WINDOW_OBSTRUCTION.md](gadget_insertion/WINDOW_OBSTRUCTION.md) |
| Window correction by a cancelling potential | Proved conditional reduction; no successful potential claimed | [CYCLE_WINDOWS.md, Section 7](gadget_insertion/CYCLE_WINDOWS.md#7-a-conditional-certificate-with-a-telescoping-correction) |
| Natural-tree paths through P128 and odd windows through order 9 | Independent numerical holdouts only | [long_path_study.json](gadget_insertion/long_path_study.json) |
| Bounded-gadget ordering and finite distance cones | Analytic reduction and checked implementation | [PROOF.md](gadget_insertion/PROOF.md), [gadget README](gadget_insertion/README.md) |

## 1. What improved in the rigorous argument

The previous first proof had the P4 total-cost deficit `10^-54` and converted
it through a mixed P2/P4 packing with factor `1/7`. The new work improves
both steps. The much larger certified local region permits sharper stability
estimates, yielding `eta=3.4e-19`; the new decomposition converts the full
P4 coefficient with factor `1/4`.

For a 3-edge-connected cubic graph, contract any perfect matching. The
resulting loopless 4-regular multigraph has every nontrivial edge cut of
size at least four. Tutte--Berge supplies a matching missing at most one
vertex; lifting it gives a partition into P4 paths and at most one P2.
A separate balanced-cut/small-shore argument transfers the same asymptotic
coefficient to every connected bridgeless cubic graph. Its error is a
constant for each fixed coefficient slack, not a sum over all decomposition
pieces. An analogous 3-edge-cut argument reduces future work further to
cyclically 4-edge-connected graphs. It does not imply long cycle lengths.

For the analytic step, nearly optimal old pairs force a P4 mass vector near
the natural tree vector. The new interval certificate controls every mass
in the full simplex with L1 distance at most `9/2500=0.0036`, including
previously absent types. Its image satisfies

\[
C_{ij}\ge (C_T)_{ij}-\frac1{300}\qquad(i\ne j).
\]

Submodular Gaussian comparison and certified integration at the lower
comparison matrix prove a total-cost gap greater than `1/10000` there.
Explicit stability estimates force any putative point closer than `eta`
to `2*M_tree` into this neighborhood, completing the contradiction.

The wide neighborhood proves the weaker bound `Phi_4<2*M_tree-1/10000`.
It should not be confused with the older, much smaller neighborhoods on
which the stronger per-vertex value `0.13614` was certified.

The combined exact calculation proves that the circuit coefficient exceeds
the original exact tree-method coefficient by more than `2.2e-18`.
[improved_bound.json](gadget_insertion/improved_bound.json) records rational
endpoints and hashes of the proof inputs. The large upstream certificates
remain dependencies and were not replayed in this work.

## 2. Updated frontier and a rigorous reason to change methods

The public [CircuitBreaker page](https://www.yukon.org/circuitbreaker),
retrieved on 2026-11-17, reports

\[
c_{\rm reported}=\frac{6341}{1341}
=4.728560775540641\ldots.
\]

It marks results provisional, pending full mathematical review. This is a
comparison with the reported exact expression, not an independent verification
of that submission. Its proof dossier was not recovered. The earlier cached
coefficient `4.672351649483386` is superseded for this comparison.
Under the same conversion, matching this coefficient requires

\[
\kappa=\frac{1341}{10000}=0.1341;
\]

beating it requires a strictly smaller `kappa`.

The certified feasible nonuniform P4 point has per-vertex cost at least

\[
L=\frac{4904795497312673}{36028797018963968}
>0.136135422304858.
\]

Therefore no uniform upper bound over this entire P4 cone can produce a
circuit coefficient exceeding

\[
1+\frac1{2L}
=\frac{22919194006794657}{4904795497312673}
=4.672813376083066\ldots.
\]

This is a ceiling on this particular proof route, not on circuit complexity
or on other graph orderings. Even the natural P4 point has a certified lower
endpoint above `0.1341`, so the obstruction is not solely an artifact of the
nonuniform candidate. The former target `0.13614` remains an intermediate
P4 question and is no longer a record target.

## 3. Long paths, window bounds, and a second obstruction

With one fixed boundary bit `b_i` per vertex, joint optimization of the
path or cycle edges has an exact combinatorial description. Put a change
edge wherever `b_i != b_(i+1)`. The split vertices form a minimum vertex
cover of this change-edge graph. Since its components are paths or an even
cycle, its cost equals the maximum matching size.

On a path the greedy matching indicators obey

\[
z_i=D_i(1-z_{i-1}),\qquad D_i=\mathbf1[b_i\ne b_{i+1}].
\]

Unrolling this recurrence through any odd number `K` gives a pointwise upper
bound involving only `K+1` consecutive bits:

\[
H_K(b_0,\ldots,b_K)
=\sum_{\ell=1}^{K}(-1)^{\ell+1}
  \mathbf1[b_{K-\ell},\ldots,b_K\text{ alternate}].
\]

For a cycle, the optimized cost is at most the sum of these window costs,
including wraparound. On bounded paths an explicit endpoint charge suffices.
This changes a potentially enormous path sign table into a short-window
Gaussian objective. The Boolean statement is proved and checked; the required
uniform Gaussian inequality is not proved. In particular, `H_K` is generally
not submodular, so the P4 covariance comparison cannot be transferred to it
without another argument. Short cycles and repeated score coordinates must
also be treated correctly.

A fresh natural-tree study used eight independent scrambled Sobol sequences
of `2^20` Gaussian samples each, exact linear-time Boolean optimization,
and an old-P2 cost control with known expectation. Its estimates are:

| Path | Natural-tree cost per vertex | Heuristic standard error |
|---|---:|---:|
| P16 | 0.134264057 | 0.000004980 |
| P24 | 0.134055498 | 0.000003146 |
| P32 | 0.133955327 | 0.000001925 |
| P48 | 0.133853162 | 0.000002066 |
| P64 | 0.133803196 | 0.000002310 |
| P128 | 0.133724384 | 0.000000823 |

For stationary matching-edge scores, the same holdouts give:

| Window | Score coordinates | Natural-tree expectation | Heuristic standard error |
|---|---:|---:|---:|
| H3 | 4 | 0.137694477 | 0.000001935 |
| H5 | 6 | 0.134044743 | 0.000001419 |
| H7 | 8 | 0.133686826 | 0.000001668 |
| H9 | 10 | 0.133652458 | 0.000002163 |

The H7 natural-tree estimate lies about `0.00041317` below `0.1341`.
That is numerical headroom at one covariance, not a universal margin;
the adversarial check below shows why it cannot justify a uniform bound.
Standard errors measure variation across independent scrambles; they are
not interval error bounds or guaranteed confidence intervals. All raw
replicate values, control means, seeds, and factorization diagnostics are
saved. An independent audit checked the marginals, controls, window indexing,
Boolean costs, and reported means and standard errors.

The conditional graph reduction is now complete: a uniform centered H7
expectation at most `beta`, together with total-cost bounds `r*beta` for
cycles of lengths `r=4,5,6,7`, implies `kappa<=beta`. For fixed large `L`,
keep cycles shorter than `L` whole and partition longer cycles into paths
of length `L` through `2L-1`. The exact path endpoint overhead is at most
5 per path. Taking `L` large before graph order removes this overhead while
retaining the bounded-gadget concentration and tie estimates. A separate
centering lemma applies to the full submodular gadget cost, so H7 itself
need not be submodular or maximized at threshold zero.

**The bare H7 hypothesis is false on the projected P8 tree face.**
[WINDOW_ADVERSARIAL.md](gadget_insertion/WINDOW_ADVERSARIAL.md) records an
independent numerical search and holdouts, including reversal averaging.
[WINDOW_OBSTRUCTION.md](gadget_insertion/WINDOW_OBSTRUCTION.md) then gives
the exact positive mass vector

\[
m=\frac1{1140}(128,640,64,160,80,40,20,2,5,1)
\]

on the ten branch-ordered P8 tree types. At that covariance, both H7
orientations and their average have expectation greater than
`17/125=0.136>0.1341`. The proof needs no eight-dimensional integration:
`H7>=D_6(1-D_5)` and its reversal reduce the lower bound to three-score
events, whose expectations are certified by the pair-sign arcsine formula.
Each event probability exceeds `0.136076509488864`.

The witness concerns the relaxed domain. It does not prove that a finite
graph realizes those masses or that an entire cycle has a bad average.
It rules out proving the desired bare window bound on this whole relaxation.

The new potential interface retains this distinction. For a fixed function
`Q` on seven ordered bits, replace H7 by

\[
W_Q(b_0,\ldots,b_7)=H_7(b_0,\ldots,b_7)
+Q(b_0,\ldots,b_6)-Q(b_1,\ldots,b_7).
\]

The corrections sum to zero around every cycle. On a path, they add at most
`max Q-min Q` to the endpoint charge. Thus a uniform centered bound
`E W_Q<=67/500`, together with the short-cycle bounds, would give
`kappa<=67/500` and `c>=317/67`, exceeding the reported `6341/1341`.
No such potential or uniform Gaussian bound is supplied here. This is a
precise next certificate interface; additional graph-realizability or
averaging constraints are another way to change the failed proof obligation.

## 4. Replay and earlier work

From the repository root:

```sh
python -B research/gadget_insertion/verify_stability_refinement.py \
  --output research/gadget_insertion/stability_refinement.json
python -B research/gadget_insertion/global_reduction_checks.py \
  --output research/gadget_insertion/global_reduction_checks.json
python -B research/gadget_insertion/global_reduction_cyclic3_checks.py \
  --output research/gadget_insertion/global_reduction_cyclic3_checks.json
python -B research/gadget_insertion/certification/mass_neighborhood.py \
  --output research/gadget_insertion/certification/mass_neighborhood.json
python -B research/gadget_insertion/certify_improved_bound.py \
  --output research/gadget_insertion/improved_bound.json
python -B research/gadget_insertion/frontier_obstruction.py \
  --output research/gadget_insertion/frontier_obstruction.json
python -B research/gadget_insertion/cycle_windows.py \
  --output research/gadget_insertion/cycle_window_checks.json
python -B research/gadget_insertion/certification/window_obstruction.py \
  --output research/gadget_insertion/certification/window_obstruction.json
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/long_path_study.py \
  --sizes 4 8 16 24 32 48 64 128 --power 20 --repeats 8 --seed 88123 \
  --output research/gadget_insertion/long_path_study.json
```

The [local certificate README](gadget_insertion/certification/README.md)
gives the wider comparison-matrix integration replay. The exact integer
[tree constant verifier](mtree_refinement/certify.py) independently encloses
`M_tree`; its derivation is in [mtree_refinement/README.md](mtree_refinement/README.md).
The [gadget README](gadget_insertion/README.md) retains P4 adversarial searches,
complete P5--P7 enumeration, finite graph inclusion checks, and the older
local certificates. These remain useful tools, with their numerical and
rigorous scopes kept explicit.

The first qualitative and quantitative proofs remain in
[GLOBAL_STRATEGY.md](gadget_insertion/GLOBAL_STRATEGY.md) and
[QUANTITATIVE_GLOBAL.md](gadget_insertion/QUANTITATIVE_GLOBAL.md).
The current result is the stronger theorem at the top of this page.
