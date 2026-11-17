# Bounded gadget insertion: reproducible research tools

This directory explores an extension of the upstream matched-edge insertion
rule to bounded paths. The exact Boolean optimizer and the tree covariance
formula are implemented separately from numerical Gaussian integration.
`PROOF.md` records the graph-theoretic argument and the remaining obligations.

**No file here certifies a better universal cubic cutwidth bound.** A low
expected cost on the infinite cubic tree does not bound every cubic graph.
Do not assume the tree covariance maximizes the relaxed distance-type cone
for a larger gadget: the numerical study below finds a nearby larger value.
All reported transformed circuit coefficients are
conditional illustrations, not established lower bounds or leaderboard
submissions.

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
are numerical evidence; they are not a rigorous comparison or an upper bound
on the cone supremum. In particular, the natural tree value is not a safe
conjectural upper bound without further work.

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

## What would be needed for a new bound

`PROOF.md` gives the combinatorial insertion extension and a concrete
graph-wide packing: choose `L=4` in a 2-factor of a triangle-free cubic kernel,
keep the perfect matching as base edges, and partition each cycle into paths
with 4 through 7 vertices. No separate short-cycle gadgets are needed for
this choice. Matching chords remain base edges; they are not freed by the
path-only Boolean optimizer.

The remaining numerical proof obligation is a **uniform, rigorous Gaussian
bound for P4, P5, P6 and P7 over every relevant neighborhood**, including
coinciding boundary keys and extra edges. This is separate from evaluating
tree points or locating a stationary point. The present deterministic
gradient code stops at P4, and explicit type enumeration is limited to P6;
these tools therefore need further development for the complete L=4 target.

The upstream theorem that the two-vertex tree point maximizes its 123-type
cone cannot be reused as a theorem about these larger gadgets. A proposed
universal constant needs a new inequality and a rigorous certificate covering
the full domain before it implies a circuit lower bound.
