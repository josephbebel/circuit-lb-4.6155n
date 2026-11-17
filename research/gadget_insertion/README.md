# Bounded gadget insertion: proofs, certificates, and further research

This directory extends the upstream matched-edge insertion rule to bounded
paths and cycles. The exact Boolean optimizer, graph reductions, rigorous
Gaussian certificates, and exploratory searches have separate scopes.
[PROOF.md](PROOF.md) gives the common ordering mechanism.

## Current rigorous result

The latest [stability refinement](STABILITY_REFINEMENT.md) and
[P4 graph reduction](GLOBAL_REDUCTION.md) prove

```
eta   = 17/(5*10^19) = 3.4e-19
M_4  <= 2*M_tree - eta
kappa <= M_4/4 <= M_tree/2 - 8.5e-20
c     >= 1 + 1/(M_tree - eta/2).
```

The exact numerical consequences, with outward rounding, are in
[improved_bound.json](improved_bound.json). In particular,

```
c >= 4.615543275582792607402913476760178553046176686.
```

This is a strict universal improvement with a very small absolute gain.
It uses the named upstream Gaussian theorem and finite certificates; the
large upstream certificates were not replayed. It does not establish a
new leaderboard record.

The graph argument first produces a P4 factor with at most one P2 in
3-edge-connected cubic graphs. A separate decomposition transfers the same
asymptotic coefficient through all 2-edge cuts, without accumulating a
linear assembly loss. Uniform P5--P7 certificates are unnecessary for this
P4 route. A further [three-edge-cut reduction](GLOBAL_REDUCTION_CYCLIC3.md)
transfers any coefficient from cyclically 4-edge-connected cubic graphs;
it does not force arbitrarily long cycles.

The original qualitative proof remains in [GLOBAL_STRATEGY.md](GLOBAL_STRATEGY.md),
and the historical `10^-54` quantitative proof remains in
[QUANTITATIVE_GLOBAL.md](QUANTITATIVE_GLOBAL.md). Their improved replacement
is the theorem above.

```sh
python -B research/gadget_insertion/verify_stability_refinement.py \
  --output research/gadget_insertion/stability_refinement.json
python -B research/gadget_insertion/global_reduction_checks.py \
  --output research/gadget_insertion/global_reduction_checks.json
python -B research/gadget_insertion/global_reduction_cyclic3_checks.py \
  --output research/gadget_insertion/global_reduction_cyclic3_checks.json
python -B research/gadget_insertion/certify_improved_bound.py \
  --output research/gadget_insertion/improved_bound.json
```

## Updated frontier and the limit of the P4 route

On 2026-11-17 the public CircuitBreaker page reported the provisional
coefficient `6341/1341=4.728560775540641...`, pending full mathematical
review. Matching it through the same circuit reduction requires
`kappa=1341/10000=0.1341`. Its underlying proof was not independently
verified here.

[FRONTIER_OBSTRUCTION.md](FRONTIER_OBSTRUCTION.md) proves that the current
uniform P4-cone route cannot reach that coefficient. One certified feasible
point already costs more than `0.136135422304858` per vertex, placing a
ceiling below `4.672813376083068` on the circuit coefficient obtainable
from a uniform bound for this particular objective. This is not a ceiling
on other methods or on circuit complexity itself. The old `0.13614` target
is now only an intermediate P4 bound.

```sh
python -B research/gadget_insertion/frontier_obstruction.py \
  --output research/gadget_insertion/frontier_obstruction.json
```

## Larger certified P4 regions

The [local certification package](certification/README.md) retains the
rigorous point enclosures:

| Exact P4 point | Certified cost per vertex, rounded outward |
|---|---:|
| Natural tree | `[0.136126085221635, 0.136126385578789]` |
| Rationalized nonuniform candidate | `[0.136135422304858, 0.136135429221129]` |

The candidate exceeds the natural point by strictly more than `9e-6`
per vertex, disproving natural-tree extremality in this relaxation.
The older entrywise covariance boxes of radius `2^-20`, and full-simplex
mass balls of L1 radius `2^-30`, satisfy the stronger local per-vertex bound
`0.13614`.

The new [wide comparison certificate](certification/WIDE_NEIGHBORHOOD.md)
proves the weaker total-cost bound `Phi_4<2*M_tree-1/10000` whenever a valid
correlation matrix satisfies `C_ij >= (C_T)_ij-1/300`. The new
[mass derivative certificate](certification/MASS_NEIGHBORHOOD.md) puts the
entire L1 mass ball of radius `9/2500=0.0036` into that comparison region.
This much wider neighborhood is what improves the global stability proof.

## Longer paths and finite windows

[CYCLE_WINDOWS.md](CYCLE_WINDOWS.md) proves an exact change-edge matching
formula for joint path/cycle optimization and pointwise odd-window upper
bounds. It also states the uniform Gaussian hypotheses needed to turn
those bounds into a universal coefficient. The hypotheses remain open.
The individual window observables are generally not submodular, so their
Gaussian bounds cannot be obtained by reusing the P4 comparison rule.

[LONG_PATH_STUDY.md](LONG_PATH_STUDY.md),
[long_path_study.py](long_path_study.py), and
[long_path_study.json](long_path_study.json) record a fresh study of natural
tree paths with 16, 24, 32, 48, 64, and 128 vertices, together with odd
stationary windows of orders 3, 5, 7, and 9. Eight independent holdouts of
`2^20` points give P32 cost about `0.133955327` and H7 cost about `0.133686826`.
These estimates lie below `0.1341` at the natural point; they do not bound
an entire cone or establish a record. Means, heuristic scramble standard
errors, raw replicates, control variates, and seeds are saved and were
independently audited.

The adversarial check is decisive for the bare window target.
[WINDOW_ADVERSARIAL.md](WINDOW_ADVERSARIAL.md) reports above-target tree-face
points, and [WINDOW_OBSTRUCTION.md](WINDOW_OBSTRUCTION.md) rigorously proves
that forward H7, reverse H7, and their average each exceed `17/125=0.136`
at one positive rational point in that relaxation. A simple three-score
lower bound and interval arcsines suffice. This is not a graph-realizability
or cycle-average claim. Section 7 of `CYCLE_WINDOWS.md` gives a conditional
alternative using a potential whose contributions cancel between neighboring
windows; no successful potential is claimed.

```sh
python -B research/gadget_insertion/cycle_windows.py \
  --output research/gadget_insertion/cycle_window_checks.json
python -B research/gadget_insertion/certification/window_obstruction.py \
  --output research/gadget_insertion/certification/window_obstruction.json
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/long_path_study.py \
  --sizes 4 8 16 24 32 48 64 128 --power 20 --repeats 8 --seed 88123 \
  --output research/gadget_insertion/long_path_study.json
```

The rest of this page documents the earlier enumerators, optimizers, and
checks, which remain available for further work.

## Coordinates and objective

A path has vertices `0,...,k-1`, internal edges `(j,j+1)`, and `k+2` boundary
incidences whose owner vertices are

```
[0, 0, 1, ..., k-2, k-1, k-1].
```

For boundary bits `s`, the objective is the minimum, over all internal edge
bits, of the number of vertices whose three incident bits do not all agree.
The routine `cost(signs)` evaluates this integer objective using a two-state
dynamic program. `least_minimizer(signs)` returns the coordinatewise least
optimal internal assignment. Both accept an array of shape `(samples,k+2)`.

`verify_boolean.py` independently enumerates every internal assignment at
each boundary pattern for small paths. It checks the dynamic-programming
values, the least minimizers, complement symmetry, all two-coordinate
submodularity inequalities, and monotonicity along every Boolean-lattice
cover relation. These finite checks support the implementation; the general
lattice argument is given in `PROOF.md`.

## Reproduce the checked finite calculations

From the repository root, using Python 3, NumPy and SciPy:

```bash
python -B research/gadget_insertion/verify_boolean.py --max-size 7
```

The checked paths have 2 through 7 vertices. The covariance check additionally
verifies the exact two-vertex tree entries and the common-process marginals
for paths through 12 vertices. Boolean arithmetic is integer arithmetic.
The covariance marginal comparison and eigendecomposition are floating-point
checks, not interval proofs.

## Reproduce the exploratory tree study

```bash
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/run_tree_study.py \
  --sizes 2 3 4 6 8 12 --power 20 --repeats 8 --seed 1729 \
  --output research/gadget_insertion/tree_study.json \
  --exact-covariances research/gadget_insertion/tree_covariances_exact.json
```

The stored run uses 8 independent scrambled Sobol sequences of 1,048,576
points each. Every path is embedded as a prefix of one common tree process,
so paths within a scramble use common random numbers. The positive leading
Gaussian factor is integrated by sorting its boundary-crossing thresholds
and summing the resulting normal-CDF interval masses. The remaining Gaussian
integral is numerical.

| Path vertices | Estimated tree cost per vertex | Heuristic standard error across scrambles |
|---:|---:|---:|
| 2 | 0.138291012 | 0.000001894 |
| 3 | 0.136923589 | 0.000001489 |
| 4 | 0.136125783 | 0.000001943 |
| 6 | 0.135300401 | 0.000002972 |
| 8 | 0.134885788 | 0.000002399 |
| 12 | 0.134471569 | 0.000001559 |

The independently known numerical two-vertex reference is
`0.1382918034411867355...`. Its difference from the estimated mean is about
`7.92e-7`. Agreement at this scale is a numerical consistency check.
The displayed standard errors measure variation across scrambles; they are
not rigorous integration errors or guaranteed confidence bounds.

The exact tree field covariance at vertex distance `d` is

```
r_d = (1+d/3) * 2^(-d/2).
```

For two distinct boundary edges with owners at path distance `d`, their
normalized key covariance is

```
(r_d + 2*r_(d+1) + r_(d+2)) / (2 + 2*r_1).
```

`exact_tree_covariance(k)` returns each entry as an exact rational number
plus an exact rational multiple of `sqrt(2)`. The JSON file stores those
rational coefficients as strings. The numerical integrator converts them
to ordinary floating point.

## API for adversarial covariance searches

```python
from path_gadgets import ConditionalGaussianCost, expected_cost, tree_covariance

# normals.shape == (number_of_samples, k+1)
estimate = expected_cost(tree_covariance(k), normals)
evaluator = ConditionalGaussianCost(tree_covariance(k))
per_sample_conditional_costs = evaluator.sample_values(normals)
```

The returned estimate is the **total cost of the gadget**, so divide by `k`
to obtain cost per vertex. Reusing standard-normal inputs across candidate
covariances can reduce noise in comparisons. A separate independent holdout
is still needed to detect optimization of numerical noise. This interface
requires a strictly positive leading factor, as supplied by the positive
distance-type correlation matrices under study.

The eigensolver allows negative eigenvalues within a small floating-point
tolerance and clips them to zero. This convention is solely for exploration.
A final certificate must establish positive semidefiniteness and integration
errors rigorously.

## Distance-type cones and an adversarial stationary-point search

`path_cone.py` enumerates relative-distance patterns around a path. It imposes
the edge Lipschitz conditions and the shortest-path descent condition at each
core vertex, treating a source at a core vertex separately. Boundary vertices
are labeled occurrences; allowing them to be independent enlarges the
relaxation and includes neighborhoods with repeated external neighbors and
extra edges. Actual coincidences remain representable by equal rows or equal
Gaussian coordinates. The cone calculation applies both row normalization
and edge-key normalization.

The enumerations contain 123 types for P2, 763 for P3, and 4,655 for P4.
`verify_cone.py` checks consistency with the upstream P2 cone and finite graph
distance profiles. Run:

```bash
python -B research/gadget_insertion/path_cone.py --sizes 2 3 4
python -B research/gadget_insertion/verify_cone.py \
  --output research/gadget_insertion/cone_verification.json
```

`path_gradient.py` computes the mixed Boolean differences' Walsh coefficients
with exact integers, then evaluates Gaussian derivatives using ordinary
Gauss–Legendre quadrature. Its P4 derivative uses a one-dimensional integral
for a conditional fourth sign moment. `cone_directional_derivatives` includes
the derivatives of both covariance normalizations. The implementation
supports P2, P3 and P4; it supplies numerical diagnostics, not error bounds.

`search_cone.py` solves numerical stationarity equations on the face supported
by the tree's outward branches. It evaluates the resulting point at two
quadrature orders, compares it with independent conditional Sobol estimates,
and checks all generator directions numerically:

```bash
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/search_cone.py \
  --size 4 --order 64 --power 18 --repeats 4 \
  --output research/gadget_insertion/cone_search_P4.json
```

The saved P4 run reports the following values at quadrature order 128:

| Cone point | Numerical expected cost per vertex |
|---|---:|
| Natural tree distribution | 0.1361262327372273 |
| Nearby stationary point on the tree face | 0.13613542578979634 |

The difference is about `9.19305e-6` per vertex. Orders 64 and 128 agree to
roughly floating-point precision in this run, and all four independent
Sobol scrambles also give the candidate a larger value. These observations
are numerical diagnostics. The independent interval certificate described
above now proves a gain greater than `9e-6`; neither calculation proves that
the candidate globally maximizes the cone.

The largest absolute support-direction derivative at the candidate is about
`5.57e-14`. The largest derivative toward an off-face generator is about
`-0.00103864`. A stationary point can fail to be a local or global maximum,
and negative outward first derivatives do not control distant points or the
entire nonlinear cone. A point in the relaxation also need not be realized
by a cubic graph.

## Independent derivative checks

```bash
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/verify_gradient.py \
  --output research/gadget_insertion/gradient_verification.json
```

The saved verification passes 344 exact Walsh inversion identities. It checks
P2 values against upstream `gc4_lib`, P2 derivatives against `bk7_core` and
`phifun`, and fourth sign moments against the upstream all-six-correlation
Plackett path. It also checks selected P4 derivatives by finite differences
of scalar point quadrature, both covariance normalizations by independent
finite differences, mass homogeneity and scale invariance, path-reversal
symmetry, and the saved candidate's support residual at order 128.

The largest observed P4 scalar finite-difference discrepancy is below
`4.45e-9`; the normalization-chain discrepancies are below `7.58e-10`.
These are implementation consistency checks at explicitly chosen points.
They are not certified integration errors, bounds over a region, or a proof
of optimality.

## P4 searches outside the initial tree face

`global_search_P4.py` uses numerical optimization over changing supports and
column generation with an adjoint through both covariance normalizations.
Every column round scans all 4,655 generators. A natural-tree mass of `10^-6`
is retained for positive definiteness, so the actual tested direction targets
are `(1-epsilon)*e_i+epsilon*m_tree`, not the unregularized vertices. All
integration and derivative calculations in this search are floating point.

```sh
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/global_search_P4.py \
  --output research/gadget_insertion/global_search_P4.json
```

All 12 saved starts returned to the six-type face. The maximum value checked
at order 128 was `0.13613542578958404` per vertex. This is evidence from a
finite search and does not cover singular points, distant unexplored
configurations, or the whole cone.

`saddle_diagnostic.py` separately evaluates a numerical Hessian at the P4
candidate in all 15 correlation coordinates. It finds four positive
eigenvalues. This diagnostic obstructs a simple ambient-concavity argument;
it does not disprove constrained optimality in the distance-type cone.

## Larger paths: efficient enumeration and adversarial experiments

The cone enumerator now supports P7. It first enumerates the core path's
`3^(k-1)` increment sequences, then forms only boundary choices satisfying the
descent predicates. For a generic type, a vertex already having a downward
path neighbor places no additional descent restriction on its boundary
leaves; otherwise at least one leaf must be one level lower. A special core
source must be the unique core minimum, has all neighbors one level higher,
and imposes the original distinctness bounds on the remaining leaves.
These products are disjoint, so enumeration needs no large deduplication set.

| Path | Generic types | Special types | Total types |
|---:|---:|---:|---:|
| P5 | 24,883 | 3,404 | 28,287 |
| P6 | 149,299 | 21,900 | 171,199 |
| P7 | 895,795 | 137,716 | 1,033,511 |

`enumerate_types(k)` retains the original sorted tuple API. For larger
calculations, `compact_types(k)` stores sorted signed-byte rows, `iter_types(k)`
streams types, and `count_types(k)` counts the exact Cartesian products.
The full compact P7 array occupies 16,536,176 bytes. `tree_support(k)` returns
only the tree's `k+2` types and exact rational masses.

```bash
python -B research/gadget_insertion/path_cone.py --sizes 5 6 7
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/verify_large_paths.py \
  --output research/gadget_insertion/large_path_verification.json
```

The saved verification establishes exact equality with the original exhaustive
enumerator for P2 through P5, checks every P7 row against the original integer
predicates, and finds no duplicate rows. It also includes all 70,816 source
distance patterns from oriented P7 paths and boundary labelings in the cube,
Petersen graph, and 7-prism. Independent covariance comparisons include repeated
neighbors and matching chords with identical boundary keys. Numerical checks
compare the new conditional-gradient evaluator with the deterministic P4
evaluator and compare its adjoint normalization calculus with the earlier
forward derivative.

`large_path_numerics.py` supplies conditional Gaussian integration for larger
Boolean tables and an efficient adjoint derivative. For optimization, it uses
a smoothly varying Cholesky representation rather than changing eigenbases:
if `alpha_star = 1/(1^T C^-1 1)`, it integrates a common Gaussian with variance
`0.95*alpha_star` and represents the positive-definite residual by Cholesky.
The derivative toward a generator can then be computed as a quadratic form,
without storing a million matrix-valued derivatives.

Reproduce the saved searches with:

```bash
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/search_large_paths.py \
  --size 5 --seed 37006 --output research/gadget_insertion/large_path_search_P5.json
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/search_large_paths.py \
  --size 6 --seed 37007 --output research/gadget_insertion/large_path_search_P6.json
OPENBLAS_NUM_THREADS=1 python -B research/gadget_insertion/search_large_paths.py \
  --size 7 --seed 37008 --output research/gadget_insertion/large_path_search_P7.json
```

Each run uses a symmetric tree-face search followed by three unconstrained
tree-face starts; estimates a gradient and scans all cone generators; tests
144 single-generator mixtures and 2,000 random sparse mixtures; locally
re-optimizes four shortlisted points; and evaluates five selected points on
six independent scrambles of `2^19` holdout points each. Training and holdout
seeds are distinct and stored in the result files.

| Path | Best tested holdout mean per vertex | Heuristic standard error | Paired gain over the natural tree |
|---:|---:|---:|---:|
| P5 | 0.135640918947 | 0.000001160 | 0.000008498 |
| P6 | 0.135312151726 | 0.000002097 | 0.000007400 |
| P7 | 0.135073733611 | 0.000001683 | 0.000006323 |

In these runs the best tested points lie on the tree face. The sparse-mixture
searches did not find a larger value. The paired gains have much smaller
heuristic standard errors, at most `7.3e-8`, because each comparison uses common
holdout normals. Neither these gains nor the displayed standard errors are
rigorous error bounds.

All tested means lie below both exploratory targets `0.13614` and `0.13615`.
This is evidence for where further certification effort may be useful:
P4 appears tighter in these experiments. It is **not** a certificate for
P5–P7. A numerical first-derivative scan at one point, finite random searches,
and local optimization do not bound the entire cone or prove that another
configuration cannot exceed either target.

## What remains for a stronger uniform bound

`PROOF.md` gives the combinatorial insertion extension and a concrete
graph-wide packing: choose `L=4` in a 2-factor of a triangle-free cubic kernel,
keep the perfect matching as base edges, and partition each cycle into paths
with 4 through 7 vertices. No separate short-cycle gadgets are needed for
this choice. Matching chords remain base edges; they are not freed by the
path-only Boolean optimizer.

That older L=4 route would require uniform rigorous Gaussian bounds for
P4, P5, P6 and P7, including coinciding boundary keys and extra edges. The
new graph decomposition removes the P5--P7 requirement for a P4-based
coefficient. However, the certified obstruction shows that even an optimal
uniform P4 bound cannot reach the updated reported frontier.

The finite-window route has a complete conditional graph reduction in
`CYCLE_WINDOWS.md`, including a version with a fixed cancelling potential.
A centered uniform bound of `67/500` for the corrected window, together
with the same per-vertex bound for cycles of lengths 4 through 7, would give
`kappa<=67/500` and `c>=317/67`. Those Gaussian inputs and a suitable
potential are unproved. The bare H7 target, even after averaging orientations,
is rigorously false on the relaxed tree face; natural-point estimates do
not establish its feasibility on that domain.

The upstream theorem that the two-vertex tree point maximizes its 123-type
cone cannot be reused as a tree-extremality theorem for larger gadgets.
The new equality/stability argument gives one weak global inequality, while
the point comparison shows natural-tree extremality itself is false for P4.
A record improvement requires a changed objective with a uniform certificate,
or additional averaging/realizability constraints that avoid the proved
obstruction to the present objective.
