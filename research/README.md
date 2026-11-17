# Cubic cutwidth research: certified refinement and larger insertion gadgets

This work starts from upstream commit
[`a56500014ba159397903dea81d41ac836f6880d5`](https://github.com/rrwilliams0/circuit-lb-4.6155n/tree/a56500014ba159397903dea81d41ac836f6880d5).
All additions are in `research/`; the upstream proof and certificates retain
their original constants.

## Results and status

| Result | Status | Where to inspect it |
|---|---|---|
| A 30-place rational enclosure of the existing tree constant | Certified by integer interval arithmetic with a proved infinite-tail bound | [mtree_refinement/README.md](mtree_refinement/README.md), [certify.py](mtree_refinement/certify.py), [certificate.json](mtree_refinement/certificate.json) |
| Insertion of bounded gadgets, with compatible monotone choices and a complete kernel cover | Analytic reduction, using the stated upstream Gaussian and circuit reductions | [gadget_insertion/PROOF.md](gadget_insertion/PROOF.md) |
| Boolean path costs and least minimizers for paths through seven vertices | Exhaustively checked against direct enumeration | [gadget_insertion/verify_boolean.py](gadget_insertion/verify_boolean.py) |
| Generalized path distance cones | Source-to-type proof plus independent finite graph and exact algebra checks | [gadget_insertion/path_cone.py](gadget_insertion/path_cone.py), [verify_cone.py](gadget_insertion/verify_cone.py) |
| Larger-path tree values and a nonuniform P4 stationary candidate | Numerical research evidence; no certified quadrature errors or global maximum | [gadget_insertion/README.md](gadget_insertion/README.md), [cone_search_P4.json](gadget_insertion/cone_search_P4.json) |

**The work does not establish a new leaderboard record or a new universal
cutwidth bound below the exact original tree constant.** It provides a
strictly stronger numerical certificate for the existing theorem and a
concrete extension whose remaining proof obligations are explicit.

## 1. Certified improvement available now

The new standalone verifier proves

\[
0.276583606882373471022763934877
\le M_{\rm tree}\le
0.276583606882373471022763934878.
\]

Combining this with the upstream theorem and reduction gives

\[
\liminf_{n\to\infty}\frac{C_{B_2}(f_n)}n
\ge 4.615543275582792605180647436566.
\]

The original certified decimal was `4.6155432392`. The new calculation uses
Python integers, rational interval endpoints, integer square roots, and
convergent series with explicit remainders. Its default run takes about a
quarter of a second in the development environment. It has no third-party
dependency. A full derivation and the arithmetic trust base are in the
refinement note.

```sh
python -B research/mtree_refinement/certify.py
```

This calculation independently certifies the explicit numerical functional.
The circuit conclusion still depends on the upstream graph-theoretic
theorems and large certificates. Those large certificates were not replayed
in this work.

## 2. Why the existing optimization cannot reach a new frontier

Upstream `proof.md`, Section 1.2, proves that its 123-type matched-insertion
supremum equals the tree value. Thus the exact circuit coefficient available
from that optimization is at most

\[
1+1/M_{\rm tree}=4.6155432755827926\ldots.
\]

Increasing floating-point precision or repeating the same cone optimization
cannot produce a substantial gain. The larger-gadget work changes the
ordering mechanism by optimizing all internal edges of a bounded path at
once, with the boundary edge scores fixed.

The cubic normalization matters:

\[
\operatorname{cw}(H)\le(\kappa+o(1))|V(H)|,
\qquad M=2\kappa,
\qquad c=1+\frac1{2\kappa}.
\]

Smaller kappa gives a larger circuit coefficient c. The original tree value
has kappa approximately `0.1382918034411867355`.

## 3. Frontier reference and useful exact targets

On 2026-11-17 the public search retrieval of
[CircuitBreaker](https://www.yukon.org/circuitbreaker) displayed the circuit
coefficient `4.672351649483386` and explicitly marked results provisional,
pending full mathematical review. The direct page reader returned a dynamic
loading shell. The leading proof, its exact rational score, author, and review
dossier were not recovered or independently verified. This snapshot can lag
subsequent submissions and should be refreshed before claiming a record.

If comparing through the same reduction, the displayed decimal corresponds
arithmetically to

\[
\kappa_{\rm display}
=\frac1{2(4.672351649483386-1)}
\approx0.136152538679224335.
\]

This is a comparison inferred from the displayed number, not a constant
verified from that submission's proof. Two exact prospective targets are:

| Universal cutwidth target, still unproved | Circuit consequence if proved |
|---|---|
| kappa <= 2723/20000 = 0.13615 | c >= 12723/2723 = 4.672420124862284... |
| kappa <= 6807/50000 = 0.13614 | c >= 31807/6807 = 4.672689878066696... |

Both would exceed that displayed decimal, subject to the full theorem and
the current official comparison contract. Neither inequality is established
by the experiments here.

## 4. The concrete remaining theorem

For a path P_k, let F_k be the minimum number of split vertices over all
internal edge-side assignments. Let C_k(m) be the boundary-score correlation
matrix of its nonnegative distance-type masses, as defined in the proof note.

A sufficient next theorem is

\[
\max_{k\in\{4,5,6,7\}}
\frac1k\sup_{m\ge0,\,\sum_t m_t=1}
\mathbb E F_k\bigl(\mathbf1[X\le0]\bigr)
\le\kappa_*,\qquad X\sim N(0,C_k(m)).
\]

The cover is already accounted for: remove a perfect matching from a
triangle-free cubic kernel and partition its remaining cycles into paths
with 4 through 7 vertices. Keep the matching and the cut cycle edges as base
edges. Coincident boundary labels and a cycle's closing edge are included;
their scores must be shared rather than independently chosen. This
specialization needs no separate short-cycle gadget.

The Boolean mechanism is justified for every size by submodularity, partial
minimization, and monotonicity of the coordinatewise least minimizers. The
proof note also supplies the truncation-shell estimate and concentration
argument. The outstanding step is the uniform numerical/analytic inequality
over each entire cone, including zero masses and singular covariances.

The current explicit enumerator handles paths through six vertices; its
default checks enumerate P2, P3, and P4 with 123, 763, and 4655 types. The
analytic definition also specifies the P7 cone, but large-cone generation
and certification need further implementation. The Boolean optimizer already
handles P7. This distinction prevents mistaking the present tools for a
completed universal certificate.

## 5. A useful obstruction discovered during the search

For P4, deterministic floating quadrature at two orders gives

| Feasible point of the relaxed cone | Numerical cost per vertex |
|---|---:|
| Original tree distribution | 0.1361262327372273 |
| Nonuniform stationary candidate on the six-type tree face | 0.1361354257897963 |

Independent conditional Sobol comparisons support the same increase. At
the candidate, numerical derivatives toward every off-face generator are
negative, and the tree-face stationarity residual is approximately `5.3e-14`.
The mass vector, all coordinates, seeds, quadrature orders, and diagnostics
are recorded in `cone_search_P4.json`.

These computations suggest that a proof asserting tree extremality for the
new relaxation would fail. They do not constitute a rigorous counterexample
without a certified error bound. Likewise, stationarity and negative outward
derivatives do not establish local or global optimality, and an abstract cone
point need not be realized by a graph.

The candidate is numerically below both prospective cutwidth targets above,
so it remains worthwhile to investigate a global certificate. A successful
proof would need either to bound a nonuniform extremizer or to use additional
realizability/averaging structure. The proof note gives a precise optional
interface for a boundary potential whose contributions cancel around the
2-factor cycles; no suitable potential is claimed here.

```sh
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/search_cone.py \
  --size 4 --output research/gadget_insertion/cone_search_P4.json
```

## 6. Verification and next work

The verified implementation checks are documented alongside their scripts:
exact Boolean tables through P7; exact agreement with the upstream 123-type
set; 57,456 actual graph source/path/label cases; 200 exact tree Gram
identities in Q(sqrt(2)); and independent numerical derivative comparisons.
They test the enumerator and search machinery and do not certify a cone
maximum. The new exact tree certificate was separately audited against its
analytic derivation and reproduced at a higher working precision.

The most valuable continuation is to turn the P4 stationary neighborhood and
its off-face directional slack into an outward-rounded certificate, while
searching the P5–P7 cones for worse configurations. Global domain coverage,
not a more favorable sampled value, is the decisive missing ingredient.

No leaderboard submission, upstream pull request, or message to another
person is part of this branch.
