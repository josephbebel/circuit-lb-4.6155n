# Reducing the universal cutwidth target to P4 alone

## Result

Let `M_4` be a proved upper bound for the total optimized P4 insertion cost
over the full distance-type cone. The bounded-gadget insertion theorem in
`PROOF.md`, together with the elementary graph reductions below, gives

\[
\operatorname{cw}(H)\le (M_4/4+o(1))|V(H)|
\]

uniformly over the simple, triangle-free, 2-connected cubic kernels required
by the upstream circuit reduction. In particular, the previously proved
inequality `M_4 <= 2 M_tree - eta` now gives

\[
\kappa\le \frac{M_{\rm tree}}2-\frac\eta4,
\qquad
c\ge1+\frac1{M_{\rm tree}-\eta/2}.
\]

The significance is structural: any proved uniform P4 bound now transfers
without uniform P5, P6, and P7 certificates. The former intermediate target
`M_4/4<=0.13614` would therefore suffice by itself for that coefficient.
It would not match the refreshed reported frontier, which corresponds to
`kappa=0.1341`. Indeed, `FRONTIER_OBSTRUCTION.md` proves that this uniform
P4-cone route cannot reach that frontier even with an optimal certificate.
This note does not infer a global bound from a stationary point or a local
certificate.

The proof has two steps. First, a 3-edge-connected cubic graph admits a
factor consisting of P4 paths and at most one P2. Second, an asymptotic
cutwidth bound on 3-edge-connected cubic graphs extends to all connected
bridgeless cubic graphs with only a constant additive loss for each fixed
slack in the leading coefficient. The second step explicitly handles
unbalanced chains of 2-edge cuts.

## 1. Matching convention and the only external matching theorem

Graphs in the assembly argument may be loopless multigraphs; degrees and
edge cuts count multiplicity. Matching size and odd-component counts are
unchanged by replacing parallel edges with a single edge, so the usual
Tutte--Berge theorem applies:

\[
\nu(Q)=\frac12\left(|V(Q)|-
 \max_{S\subseteq V(Q)}\{o(Q-S)-|S|\}\right).
\]

Here `o` counts components of odd order. A primary exposition with a proof
is Michel Goemans's MIT notes, *Lecture notes on non-bipartite matching*,
Theorem 2.1:

<https://math.mit.edu/~goemans/18433S15/matching-nonbip-notes.pdf>.

The existence of a perfect matching in a bridgeless cubic graph is the
standard Petersen theorem already used in the upstream proof. No separate
P4-factor theorem from the literature is needed below.

## 2. A P4 factor, with at most two leftover vertices, in 3-edge-connected graphs

Let `H` be a connected loopless 3-edge-connected cubic multigraph. Except
for the two-vertex triple-edge graph, it is simple: a pair of parallel edges
in a larger connected cubic graph would make their two endpoints a shore
of a 2-edge cut. The two-vertex exception has cutwidth three and can be
handled as a bounded exceptional graph.

Choose a perfect matching `M` of a simple such `H` and contract its edges.
The resulting multigraph `Q` is connected, loopless, and 4-regular. Every
cut of `Q` lifts to a cut of `H`, hence has size at least three. Every cut
of an even-regular graph has even size, so every nontrivial cut of `Q` has
size at least four.

For any nonempty `S` properly contained in `V(Q)`, each component of `Q-S`
has at least four incident cut edges. Counting those belonging to the odd
components gives

\[
4o(Q-S)\le4|S|.
\]

The same inequality is trivial for `S=V(Q)`. At `S=empty`, connectedness
gives deficiency zero when `|Q|` is even and one when it is odd.
Tutte--Berge therefore supplies a matching `J` leaving at most one vertex
of `Q` uncovered.

Each edge of `J` lifts to an edge of `H` joining two disjoint matching
edges of `M`; their union is a path on four distinct vertices. Distinct
edges of `J` give vertex-disjoint P4 paths. If a vertex of `Q` is uncovered,
the corresponding original matching edge is the single P2 component.
Additional chords of the four vertices remain base edges, exactly as in
`PROOF.md`. Thus the graph is covered by

\[
\lfloor |V(H)|/4\rfloor\text{ P4 gadgets and at most one P2 gadget}.
\]

The bounded-gadget insertion theorem consequently bounds its expected cost
by `M_4 |V(H)|/4+O(1)`, before the uniform vanishing shell and concentration
terms. The P4 theorem in `PROOF.md` applies to simple cubic graphs; it does
not require the P4 paths to be induced or the graph to be triangle-free.
Hence every fixed `alpha>M_4/4` has a constant `C_alpha` such that

\[
\operatorname{cw}(H)\le\alpha|V(H)|+C_\alpha
\tag{1}
\]

for all 3-edge-connected cubic graphs, including the bounded two-vertex
multigraph exception after increasing `C_alpha`.

## 3. Splitting a 2-edge cut preserves the class

Let `H` be a connected bridgeless loopless cubic multigraph, and let
`delta(A)={a_1b_1,a_2b_2}`, with `a_i in A` and `b_i outside A`.
Both shores are connected: two components on one shore would each need
at least two edges leaving them, contradicting the cut size two.

The endpoints on either shore are distinct. If both cut edges met one
vertex `a`, its third edge would be the only connection of `A-{a}` to the
rest of `H`, and so would be a bridge. The case `A={a}` cannot occur in a
loopless cubic graph with cut size two. The other shore is identical.

Define the completed graphs

\[
H_A=H[A]+a_1a_2,
\qquad H_B=H[V(H)\setminus A]+b_1b_2.
\]

The added edges are called virtual edges. Both completed graphs are
connected, loopless, cubic, and bridgeless. For bridgelessness, the virtual
edge cannot be a bridge because its shore was connected. A bridge among
the original edges would separate the two virtual endpoints together on
one side, leaving a component with no boundary attachment to the other
shore of `H`; that same edge would then be a bridge of `H`.

Concatenating any layouts of `H_A` and `H_B` proves

\[
\operatorname{cw}(H)
\le\max\{\operatorname{cw}(H_A),\operatorname{cw}(H_B)\}+2.
\tag{2}
\]

For a cut inside one completed shore, deleting its virtual edge and
restoring its two external edges increases the crossing count by at most
two. The cut between the two shores has size two. Repeatedly using (2)
without further organization could accumulate a linear error; the next
two sections avoid this problem.

## 4. A core when every 2-edge cut has a small shore

Fix an integer `T>=2`, and suppose `n=|V(H)|>3T` and no 2-edge cut has
both shores of size at least `T`. Let `F` be the family of shores `A` with
`|delta(A)|=2` and `|A|<T`, and take its inclusion-maximal members.

### 4.1 Maximal small shores are disjoint and nonadjacent

If `A,B in F` intersect and neither contains the other, their union has
size less than `2T`, hence is a proper subset whose complement has more
than `T` vertices. Cut submodularity gives

\[
|\delta(A\cup B)|+|\delta(A\cap B)|
\le|\delta(A)|+|\delta(B)|=4.
\]

Both cuts have size at least two by connectedness and bridgelessness.
They therefore both have size two. The union must have size less than
`T`, since otherwise it would be a forbidden 2-edge cut with two large
shores. It is consequently in `F`. Distinct inclusion-maximal members
cannot intersect.

If two distinct maximal members are adjacent, their disjoint union has
cut size `4-2e(A,B)<=2`. It is nonempty and proper, so its cut size is
exactly two. Again the union has size less than `2T`, and the absence of a
balanced cut forces it to have size less than `T`. This contradicts
maximality. Hence no edge joins two distinct maximal small shores.

Let `X` be the vertices outside all these maximal shores. It is nonempty:
otherwise the disjoint, nonadjacent shores would disconnect `H` (and a
single shore cannot contain all `n` vertices).

### 4.2 Replacement gives a 3-edge-connected cubic core

The two external neighbors of each maximal shore lie in `X`, by the
nonadjacency just proved; they are distinct by Section 3. Replace that
shore and its two attachment edges by one virtual edge joining those
neighbors. This produces a connected loopless cubic multigraph `K` on `X`.

Every nontrivial cut of `K` lifts to a cut of `H` of the same size. Include
a replaced shore when both endpoints of its virtual edge are included,
exclude it when both are excluded, and choose either side when the
virtual edge crosses. A crossing virtual edge then contributes exactly
one crossing attachment edge in `H`.

In particular `K` has no bridge. If it had a 2-edge cut, its lift would
have a shore of size less than `T`. Such a shore belongs to `F`, and hence
is contained in a maximal member of `F`. But both lifted shores contain
at least one vertex of `X`, contradicting that containment. Thus `K` is
3-edge-connected, with the possible two-vertex triple-edge exception.

### 4.3 Expanding the small shores costs only a constant

Take a layout of `K`. For each virtual edge `ab`, place all vertices of
its associated small shore immediately after the earlier of `a,b` in the
layout, in any order. If several shores have that earlier endpoint, place
them consecutively in any order.

Between expansions, every replaced virtual edge has exactly the same
crossing contribution as its two attachment edges and its entirely
processed or unprocessed shore. During an expansion, only that one shore
is partially processed. If its size is `s<T`, it has `(3s-2)/2` internal
edges. These edges and its two attachments contribute at most
`(3s-2)/2+2`, replacing one crossing virtual edge in the core layout.
Therefore

\[
\operatorname{cw}(H)\le\operatorname{cw}(K)+\frac{3T}{2}.
\tag{3}
\]

The number of small shores does not occur in this error bound.

## 5. Uniform assembly theorem

**Theorem.** Suppose that for some `alpha>0` and `C_alpha>=0`, every
connected loopless 3-edge-connected cubic graph obeys (1). Then every
connected bridgeless loopless cubic graph obeys

\[
\operatorname{cw}(H)\le\alpha|V(H)|+C_\alpha+5T,
\qquad T\ge\max(2,\lceil2/\alpha\rceil).
\tag{4}
\]

**Proof by induction on the order.** For `n<=3T`, the trivial bound
`cw(H)<=|E(H)|=3n/2<=9T/2` proves (4).

If a 2-edge cut has two shores of size at least `T`, use Section 3 and the
induction hypothesis on its completed shores. Equation (2) gives

\[
\begin{aligned}
\operatorname{cw}(H)
&\le\alpha\max(|A|,|B|)+C_\alpha+5T+2\\
&\le\alpha(n-T)+C_\alpha+5T+2\\
&\le\alpha n+C_\alpha+5T,
\end{aligned}
\]

using `alpha T>=2`.

Otherwise Section 4 supplies the 3-edge-connected core `K`. Apply (1) and
(3) to obtain

\[
\operatorname{cw}(H)
\le\alpha|V(K)|+C_\alpha+3T/2
\le\alpha n+C_\alpha+5T.
\]

This closes the induction. No bound on the depth of an arbitrary
2-edge-cut decomposition is assumed. **QED.**

Applying this theorem separately for every `alpha>M_4/4` proves the stated
uniform asymptotic coefficient `M_4/4` for all required kernels.

## 6. A sharper direct P2/P4 packing, and its limitation

For completeness, even without the assembly theorem, replacing the old
maximal matching by a maximum matching improves its packing guarantee.
Let `Q` be any connected loopless 4-regular multigraph of order `q>=7`.
For a nonempty vertex set `S`, let `a` count singleton odd components of
`Q-S`, and `b` count its other odd components. The former have four
outgoing edges. Each latter component has at least three vertices and an
even positive cut, hence at least two outgoing edges. Writing `s=|S|`,

\[
4a+2b\le4s,
\qquad q\ge s+a+3b.
\]

If `d=a+b-s`, the first inequality gives `b>=2d`, and substitution in the
second gives `q>=2a+4b-d>=7d`. For `S=empty`, the deficiency is at most
one, also at most `q/7`. Tutte--Berge yields

\[
\nu(Q)\ge3q/7.
\]

Thus a fixed perfect matching of `H` can always be extended to at least
`3|V(H)|/14` disjoint P4 gadgets, up to a bounded small-order exception.
With `delta=2M_tree-M_4`, this alone gives

\[
\kappa\le M_{\rm tree}/2-3\delta/14.
\]

The universal assembly result is stronger: it gives
`kappa<=M_tree/2-delta/4`.

The coefficient `3/7` for a fixed contracted matching is asymptotically
sharp. Start with a connected 4-regular loopless graph on `s` vertices,
and replace each of its `2s` edges by a three-vertex gadget having five
internal edges (two doubled sides and one single side of a triangle) and
one attachment at each endpoint of the single side. The result has `7s`
vertices. Deleting the original `s` vertices leaves `2s` odd components,
so every matching has size at most `3s`; the lower bound attains this.
Such examples can be lifted to simple triangle-free bridgeless cubic
graphs by orienting each doubled pair oppositely, using an Eulerian
orientation of the original graph, and splitting each vertex into an
in/out matched pair. An oriented edge `u->v` becomes an edge from the out
half of `u` to the in half of `v`, and every matched pair joins its two
halves. The result is simple and bipartite. The balanced orientation
decomposes into directed cycles; their alternating lifts cover every
edge, including every pair edge, so the lift is bridgeless. This limitation concerns the fixed-matching packing
scheme, and does not contradict the assembly theorem.

## 7. Verification scope

`global_reduction_checks.py` exercises the factor construction, 2-edge-cut
completions, maximal-small-shore core, and expansion inequality on exact
finite graphs, including multigraphs, long 2-sum chains, and the sharp
fixed-matching examples. These are implementation checks supporting the
proof; the argument in Sections 2--5 establishes the general quantifiers.
No Gaussian numerical search or unproved cone-extremality assertion enters
the structural reduction.

The stored run checks 43 graphs, 108 completed cut shores, 34
3-edge-connected factor cases, and 55 small-shore-core parameter cases.
The expanded-cube case replaces every edge of an eight-vertex cube by a
six-vertex two-terminal gadget: its 80 vertices reduce back to the cube,
and all tested expanded layouts obey (3). Long chains exercise the
balanced-cut alternative. The sharp contracted-matching examples have
21, 28, and 49 vertices and matching sizes 9, 12, and 21 respectively.

Run from the repository root with Python, NumPy, and SciPy:

```sh
python -B research/gadget_insertion/global_reduction_checks.py \
  --output research/gadget_insertion/global_reduction_checks.json
```

The MILP solver is used only to propose matching witnesses. The program
checks vertex-disjointness and the required cardinality on the exact
integer edge lists. The relevant upper bounds follow separately from
vertex count or the displayed Tutte barrier; no floating-point solver
optimality report is used as a proof input.
