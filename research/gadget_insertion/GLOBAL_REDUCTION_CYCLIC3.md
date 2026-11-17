# Removing nontrivial three-edge cuts without changing the asymptotic coefficient

This note supplements `GLOBAL_REDUCTION.md`. Its theorem is structural: it
reduces the graph class on which a positive linear cutwidth coefficient must
be proved. It supplies no new Gaussian inequality or numerical coefficient.
In particular, a cyclically 4-edge-connected cubic graph can have 4-cycles;
the theorem does not provide a factor consisting of arbitrarily long paths
or cycles.

Throughout, graphs are finite, connected, loopless cubic multigraphs. Edges
are counted with multiplicity. A graph is **3-edge-connected** if every
nonempty proper vertex set has at least three boundary edges. A 3-edge cut
is **nontrivial** if both shores have at least two vertices. In this graph
class, saying that there is no nontrivial 3-edge cut is equivalent to saying
that the graph is cyclically 4-edge-connected. Indeed, a shore of a 3-edge
cut is connected, has odd order s, and has (3s-3)/2 internal edges. When
s>1 it therefore contains a cycle; when s=1 it contains no cycle.

The graph with two vertices and three parallel edges is permitted as a
bounded exception. Every other 3-edge-connected loopless cubic multigraph
is simple: endpoints of a double edge would otherwise form a 2-edge cut.

## The coefficient-transfer theorem

Suppose that alpha>0 and C>=0 satisfy

    cw(K) <= alpha |V(K)| + C

for every cyclically 4-edge-connected loopless cubic multigraph K, including
the two-vertex exception. Choose an integer T>=2 with

    alpha (T-1) >= 3.

Then every 3-edge-connected loopless cubic multigraph H satisfies

    cw(H) <= alpha |V(H)| + C + 5T.                         (1)

Consequently any uniform asymptotic coefficient on the cyclically
4-edge-connected class transfers unchanged to the 3-edge-connected class.
Combining this with the two-edge-cut reduction in `GLOBAL_REDUCTION.md`
transfers the same coefficient to all connected bridgeless cubic graphs
needed there.

### Completing a shore

Let A be a shore of a 3-edge cut in H. Form H_A by retaining H[A], adding
one new vertex a, and replacing each of the three cut edges by an edge from
its endpoint in A to a. The graph H_A is loopless cubic and
3-edge-connected. To check the latter, take any nontrivial cut of H_A and
choose its shore S which does not contain a. Its boundary size in H_A
equals its boundary size in H, so it is at least three.

The shore A is connected: otherwise each component would have at least
three edges leaving it in H, exceeding the three available boundary edges.
If |A|>1, the three boundary endpoints in A are distinct. If two were at x,
then A\{x} would have boundary size at most two, contrary to
3-edge-connectivity. When A is a singleton, H_A is just the two-vertex
triple-edge graph.

For a 3-edge cut with complementary shores A and B, take layouts of H_A
and H_B, delete their new vertices from the layouts, and concatenate the
remaining orders. At a cut within one shore, its internal crossing edges
are a subset of those counted at a corresponding cut of that completed
layout; at most three original edges between A and B additionally cross.
At the join only those three edges cross. Thus

    cw(H) <= max(cw(H_A), cw(H_B)) + 3.                    (2)

### Induction: a balanced cut

We prove (1) by induction on n=|V(H)|. When n<=3T, the trivial bound
cw(H)<=3n/2<=9T/2 proves (1).

Suppose n>3T. If some 3-edge cut has both shores of size at least T, the
completed shores each have at most n-T+1<n vertices. The induction
hypothesis and (2) give

    cw(H) <= alpha(n-T+1) + C + 5T + 3
          <= alpha n + C + 5T.

The last step is exactly the choice alpha(T-1)>=3.

### Induction: all cuts have a small shore

It remains to consider n>3T when every 3-edge cut has a shore of size less
than T. Let S be the family of shores A with |A|<T and |delta(A)|=3, and
let A_1,...,A_k be its inclusion-maximal members.

These maximal members are pairwise disjoint. Otherwise two distinct
maximal members A,B have nonempty intersection. Their union is proper,
because |A union B|<2T<n. Cut submodularity and 3-edge-connectivity give

    6 = |delta(A)|+|delta(B)|
      >= |delta(A intersection B)|+|delta(A union B)| >= 6.

In particular, |delta(A union B)|=3. Its complement has size greater than
T, since n>3T and |A union B|<2T. The absence of a balanced 3-edge cut
therefore forces |A union B|<T. This contradicts maximality.

The members also cover V(H): every singleton is a shore with boundary
size three and belongs to S. Thus they form a partition into nonempty
sets of size less than T. Unlike the two-edge-cut construction, distinct
members may be adjacent, and no nonadjacency claim is needed.

Contract each A_i to one vertex, discard its internal edges, and retain
all edges between different parts. The resulting graph K is connected,
loopless and cubic. Every cut of K lifts, by taking the union of its
parts, to a cut of the same size in H. Thus K is 3-edge-connected.

Moreover, K has no nontrivial 3-edge cut. Such a cut would lift to a
3-edge cut of H with each shore containing at least two whole parts.
One lifted shore U has size less than T and so belongs to S. It must be
contained in an inclusion-maximal member of S, but it already contains
two distinct maximal members, which is impossible. Hence K is cyclically
4-edge-connected (or is the permitted two-vertex exception).

Take a layout of K of width at most alpha k+C. Replace its vertex i
consecutively by all vertices of A_i, in any fixed order. Between completed
parts the crossing-edge count is exactly the corresponding count in K.
At a cut inside the active part A_i of size s<T, compare with the cut of
K just before its vertex i. Crossings not incident to A_i agree. At most
three attachment edges of A_i can add crossings, and at most
(3s-3)/2 internal edges cross. Therefore the increase is at most

    3 + (3s-3)/2 = (3s+3)/2 <= 3T/2.

Only one part is active at a time, so this error is paid once per cut,
not once per part or decomposition level. Since k<=n,

    cw(H) <= cw(K)+3T/2 <= alpha n+C+3T/2,

which proves (1) and completes the induction.

## Uniform asymptotic form

If cw(K)<=(beta+o(1))|V(K)| holds uniformly on the smaller class, fix any
alpha>beta with alpha>0. Uniformity and the trivial bound on finitely many
small orders provide a constant C for the theorem. Applying (1), and then
the two-edge-cut theorem when needed, gives cw(H)<=alpha|V(H)|+O_alpha(1).
Letting alpha decrease to beta proves the claimed transfer. The constants
may depend on alpha; they never depend on H.

## Finite checks

`global_reduction_cyclic3_checks.py` enumerates all three-edge cuts of its
test graphs, verifies both shore completions and concatenation inequalities,
constructs the maximal-small-shore quotient whenever the hypotheses apply,
and verifies the expansion inequality for three layouts. Its graph families
include standard cubic examples, seeded random cubic graphs, three-edge
sums, and replacements of every cube vertex by a triangle or a seven-vertex
cube fragment. All graph and width checks use exact integer arithmetic.

Run from the repository root:

```sh
python -B research/gadget_insertion/global_reduction_cyclic3_checks.py \
  --output research/gadget_insertion/global_reduction_cyclic3_checks.json
```

The finite tests are regression checks for the constructions. The preceding
argument proves the universal theorem independently of those tests.
