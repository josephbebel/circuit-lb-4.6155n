# Cubic cutwidth research: a strict global improvement and larger-gadget certificates

This work starts from upstream commit
[`a56500014ba159397903dea81d41ac836f6880d5`](https://github.com/rrwilliams0/circuit-lb-4.6155n/tree/a56500014ba159397903dea81d41ac836f6880d5).
All additions are in `research/`; the upstream proof and certificates retain
their original constants.

## Results and status

| Result | Status | Where to inspect it |
|---|---|---|
| A 30-place rational enclosure of the existing tree constant | Certified by integer interval arithmetic with a proved infinite-tail bound | [mtree_refinement/README.md](mtree_refinement/README.md), [certify.py](mtree_refinement/certify.py), [certificate.json](mtree_refinement/certificate.json) |
| A strict universal improvement using P2/P4 insertion | Analytic proof and explicit gap `10^-54`, using the upstream finite certificates and the new local certificate | [GLOBAL_STRATEGY.md](gadget_insertion/GLOBAL_STRATEGY.md), [QUANTITATIVE_GLOBAL.md](gadget_insertion/QUANTITATIVE_GLOBAL.md), [quantitative_verification.json](gadget_insertion/quantitative_verification.json) |
| P4 point and positive-radius neighborhood bounds below `0.13614` | Certified interval integration, with a rigorous comparison disproving natural-tree extremality in the relaxation | [certification/README.md](gadget_insertion/certification/README.md), [comparison.json](gadget_insertion/certification/comparison.json) |
| Insertion of bounded gadgets, with compatible monotone choices and a complete kernel cover | Analytic reduction, using the stated upstream Gaussian and circuit reductions | [gadget_insertion/PROOF.md](gadget_insertion/PROOF.md) |
| Boolean path costs and least minimizers for paths through seven vertices | Exhaustively checked against direct enumeration | [gadget_insertion/verify_boolean.py](gadget_insertion/verify_boolean.py) |
| Generalized path distance cones | Source-to-type proof plus independent finite graph and exact algebra checks | [gadget_insertion/path_cone.py](gadget_insertion/path_cone.py), [verify_cone.py](gadget_insertion/verify_cone.py) |
| P4 adversarial search and complete P5-P7 enumeration/search | Numerical research evidence; no record-target bound over an entire cone | [gadget_insertion/README.md](gadget_insertion/README.md), [global_search_P4.json](gadget_insertion/global_search_P4.json), [large_path_verification.json](gadget_insertion/large_path_verification.json) |

**The new proof gives a universal cutwidth coefficient strictly below the
exact original tree-method coefficient. Its explicit gain is extremely
small and does not establish a new leaderboard record.** The much stronger
prospective target `kappa <= 0.13614` is certified on two local neighborhoods;
its global proof remains open.

## 1. New universal theorem

Let `M=M_tree` be the exact old pair constant and let `M_4` denote the
maximum total four-vertex insertion cost over all 4,655 distance types.
The new argument proves

\[
M_4\le 2M-10^{-54}.
\]

It first proves strictness without numerical integration. Decomposing a P4
into two old pairs gives a sum of two nonnegative pair deficits and a
nonnegative joint-insertion gain. Simultaneous equality of the old pair
bounds forces the unique natural P4 mass vector, where an explicit open
Gaussian event gives strictly positive insertion gain. Compactness excludes
equality throughout the closed relaxed cone.

The quantitative companion uses the upstream local quadratic bound and
finite cover to control a nearly optimal pair covariance. Quantitative
versions of the equality constraints then put any P4 point with two nearly
optimal pairs inside the newly certified natural-tree neighborhood. This
gives the explicit, deliberately conservative gap `10^-54`.

A maximal matching after contracting the kernel's original perfect matching
provides at least `h/7` disjoint P4 gadgets on `h` vertices; all remaining
vertices lie in P2 gadgets. The resulting universal asymptotic bound is

\[
\operatorname{cw}(H)
\le\left(\frac M2-\frac{10^{-54}}7+o(1)\right)|V(H)|.
\]

This applies to the simple, triangle-free, 2-connected cubic kernels needed
by the upstream reduction; the packing itself holds for every simple
bridgeless cubic graph.

The upstream kernel-to-circuit reduction therefore gives

\[
\liminf_{n\to\infty}\frac{C_{B_2}(f_n)}n
\ge 1+\frac1{M-2\cdot10^{-54}/7}.
\]

The exact arithmetic verifier reruns the independent tree enclosure at
120 working decimal places and proves that this circuit coefficient exceeds
the old exact method coefficient by more than `3*10^-54`. It stores complete
rational bounds; this very small improvement has no effect on the familiar
displayed digits `4.6155432755827926...`.

```sh
python -B research/gadget_insertion/verify_rigidity.py
python -B research/gadget_insertion/verify_quantitative.py \
  --output research/gadget_insertion/quantitative_verification.json
```

These exact finite checks support the written proof; they do not replace
its analytic arguments or replay the large upstream certificates.

## 2. Certified refinement of the original numerical constant

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

## 3. Why the existing optimization cannot reach a new frontier

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

## 4. Frontier reference and useful exact targets

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

## 5. The concrete remaining theorem for a competitive record

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

The factorized enumerator now handles P7 in a compact array. Its exact counts
are 4,655, 28,287, 171,199, and 1,033,511 for P4, P5, P6, and P7 respectively.
Complete generation, finite graph inclusion checks, and numerical searches
are implemented. The outstanding task is the substantially sharper global
inequality above; the tiny universal theorem in Section 1 does not imply it.

## 6. Rigorous P4 local bounds and the failure of tree extremality

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

The new interval calculations make the comparison rigorous:

| Exact cone point | Certified per-vertex enclosure, rounded outward |
|---|---:|
| Natural tree | `[0.136126085221635, 0.136126385578789]` |
| Rationalized nonuniform candidate | `[0.136135422304858, 0.136135429221129]` |

The candidate exceeds the natural tree by strictly more than `9*10^-6`
per vertex. **Natural-tree extremality is false for this relaxed P4 cone.**
An abstract cone point need not be realized by a finite graph, so this is
not a statement that a corresponding graph is extremal.

For both exact centers, every unit-diagonal symmetric covariance within
entrywise off-diagonal radius `2^-20` is positive definite and has cost per
vertex below `6807/50000=0.13614`. Each certified region also contains the
image of the full-simplex mass ball of L1 radius `2^-30`, including masses
outside the six-type face. The proof uses exact Boolean coefficients,
two-dimensional Plackett integration, outward interval Taylor remainders,
interval Cholesky, and explicit derivative bounds. See the
[certification note](gadget_insertion/certification/README.md) for replay
commands and arithmetic assumptions.

Stationarity and negative outward derivatives still do not establish local
or global optimality of the nonuniform candidate. A numerical Hessian
diagnostic also finds positive directions in the ambient correlation
coordinates, so a general concavity argument is unavailable.

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

## 7. Adversarial searches and verification

A 12-start P4 column-generation search optimized over changing supports and
checked all 4,655 regularized generator directions at each column round.
All starts returned to the six-type face, with maximum order-128 numerical
value `0.13613542578958404`. The `10^-6` natural-tree regularizer keeps those
numerical evaluations positive definite; the search does not cover singular
points or constitute a global upper certificate.

The P5-P7 searches optimized tree-face candidates, screened every generator
direction numerically, checked 2,144 sparse mixtures per size, polished four
candidates, and used six independent holdout scrambles:

| Path | Best tested holdout mean per vertex | Heuristic standard error |
|---|---:|---:|
| P5 | 0.135640918947 | 0.000001160 |
| P6 | 0.135312151726 | 0.000002097 |
| P7 | 0.135073733611 | 0.000001683 |

No tested sparse mixture exceeded the best tree-face point. These are search
outcomes, not rigorous upper bounds or guaranteed confidence intervals.

The verified implementation checks are documented alongside their scripts:
exact Boolean tables through P7; exact agreement with the upstream 123-type
set; 57,456 actual graph source/path/label cases; 200 exact tree Gram
identities in Q(sqrt(2)); and independent numerical derivative comparisons.
The larger enumeration adds exact agreement with the original brute-force
definition through P5, checks all 1,033,511 P7 rows and duplicate absence,
and includes 70,816 further actual graph source patterns. These tests verify
the enumerator and search machinery; they do not bound an entire cone.

The interval evaluator passed independent exact Taylor-coefficient and
polynomial-remainder checks, complete point replays, and a separate analytic
audit. The global proof and its constants received a separate review. The
upstream graph-theoretic results and large finite certificates remain named
dependencies and were not replayed in this work.

The decisive next mathematical task is a strong global certificate covering
the complement of the P4 local region, followed by the P5-P7 target bounds.
The new rigidity/stability argument gives one complete but very weak global
bound. Improving that argument or reducing the remaining domains is a more
useful direction than simply increasing the precision of the sampled values.

No leaderboard submission, upstream pull request, or message to another
person is part of this branch.
