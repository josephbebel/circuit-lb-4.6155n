# The refreshed frontier rules out the current uniform P4 target

## Reported frontier and exact normalization

The public [CircuitBreaker page](https://www.yukon.org/circuitbreaker),
retrieved on 2026-11-17, reports the leading linear circuit coefficient

\[
c_{\rm reported}=\frac{6341}{1341}
=4.728560775540641\ldots.
\]

The page labels its results provisional, pending full mathematical review.
This note uses the displayed exact expression for comparison; it does not
independently verify the leader's proof. A prior cached retrieval had shown
`4.672351649483386`; that snapshot is superseded for this comparison.

Under the same arithmetic conversion used by this project,

\[
c=1+\frac1{2\kappa},
\]

matching the reported coefficient requires

\[
\kappa=\frac1{2(c_{\rm reported}-1)}
=\frac{1341}{10000}=0.1341.
\]

Exceeding the reported coefficient requires `kappa<0.1341`. The former
prospective target `kappa=0.13614` would give `c=31807/6807`, which is below
the refreshed reported frontier.

## A rigorous obstruction to the present uniform-cone method

Let `M_4` denote the supremum of total P4 insertion cost over the current
4,655-type relaxed cone. The separately certified rational nonuniform point
has expected cost per vertex at least

\[
L=\frac{4904795497312673}{36028797018963968}
=0.136135422304858\ldots.
\]

Thus `M_4/4>=L>0.1341`. Any proof which establishes a uniform inequality

\[
\Phi_{P_4}(C_4(m))\le4\kappa\qquad\text{for every mass }m
\]

must have `kappa>=L`. Applying the circuit conversion to such a uniform
bound can therefore yield a coefficient no greater than

\[
\boxed{
1+\frac1{2L}
=\frac{22919194006794657}{4904795497312673}
<4.672813376083068.
}
\]

This ceiling is below `6341/1341` by more than `0.055747399457574`.
It does not depend on proving global optimality of the nonuniform point:
one rigorously evaluated feasible point suffices.

The natural P4 point also has a certified per-vertex lower endpoint above
`0.1341`. Consequently, merely replacing the nonuniform conjectural
extremizer by the natural tree value would still not meet this frontier.

`frontier_obstruction.py` checks every displayed comparison using exact
rational arithmetic. It verifies that the candidate masses are positive,
sum to one, and use types in the declared cone. It reads the point
certificate's exact lower endpoint and records the input hashes. Its output
is `frontier_obstruction.json`; it performs no additional Gaussian integration.

## What this obstruction does and does not say

This is an obstruction to the specified **uniform relaxed-P4 proof route**.
It is not an upper bound on Boolean circuit complexity, and it is not an
upper bound on what every possible graph ordering can achieve. The abstract
mass vector need not be realizable as a finite graph neighborhood.

The new `GLOBAL_REDUCTION.md` remains useful: it proves that a global P4
upper bound transfers to the kernel coefficient without P5-P7 certificates.
However, even an optimal solution of the current P4 relaxation cannot beat
the newly reported coefficient through this conversion. A new record now
requires a changed proof obligation or a changed ordering mechanism.

Concrete directions are:

1. **Optimize longer paths or whole cycles.** Joint optimization can remove
   costs retained by four-vertex gadgets. For a path with one fixed boundary
   bit at each vertex, the optimized split count equals the minimum vertex
   cover of the path edges joining unequal adjacent bits, hence the maximum
   matching size of those change edges. This leads to finite-window
   alternating-sign observables, potentially avoiding enumeration of an
   enormous long-path cone. `CYCLE_WINDOWS.md` proves the Boolean identities
   and a complete conditional bounded-gadget reduction. However,
   `WINDOW_OBSTRUCTION.md` rigorously rules out the bare H7 uniform bound
   on the projected P8 tree face, even after averaging orientations.
   A cancelling potential, justified realizability restrictions, or another
   average bound is needed to overcome that particular obstruction.
2. **Use richer bounded gadgets.** Additional internal edges or branching
   shapes change the Boolean cost functional. Their own distance cones,
   global coverage, and Gaussian certificates are then required. The P4
   point certificate places no ceiling on those different functionals.
3. **Prove stronger realizability or global averaging inequalities.** A
   point outside the realizable graph family may be removable, and a
   sum-of-gadgets bound can be stronger than an independent bound for every
   gadget. Any such argument must state and prove its additional constraints.
   Simply imposing tree extremality is inadequate: the certified natural P4
   value itself exceeds the new target.

The immediate practical implication is to stop treating `0.13614` as a
record target. The existing P4 certificate remains a proved methodological
advance and a reusable local tool, while the next competitive search needs
to address the smaller threshold `0.1341` with a different mechanism.

## Replay

```sh
python -B research/gadget_insertion/frontier_obstruction.py \
  --output research/gadget_insertion/frontier_obstruction.json
```
