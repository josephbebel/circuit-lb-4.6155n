#!/usr/bin/env python3
"""Exact finite checks used by GLOBAL_STRATEGY.md's equality analysis.

No Gaussian numerical integration is used. Distance-type membership is
integer arithmetic. The diagonal constraints, their rank, and their unique
solution are checked using fractions.Fraction.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path

from path_cone import enumerate_types, structure


def power_two(exponent):
    return Fraction(2**exponent) if exponent >= 0 else Fraction(1, 2**(-exponent))


def exact_rref(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    pivots = []
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x/scale for x in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][col]:
                factor = a[i][col]
                a[i] = [x-factor*y for x, y in zip(a[i], a[rank])]
        pivots.append(col)
        rank += 1
        if rank == len(a):
            break
    return a, pivots


def verify():
    k = 4
    types = enumerate_types(k)
    pair_types = set(enumerate_types(2))
    projections = ((0, 1, 4, 5, 2, 6), (2, 3, 1, 7, 8, 9))
    projection_results = []
    for projection in projections:
        images = [tuple(t[i]-t[projection[0]] for i in projection) for t in types]
        assert all(image in pair_types for image in images)
        projection_results.append({
            "label_indices": projection,
            "source_types_checked": len(images),
            "distinct_pair_type_images": len(set(images)),
            "all_images_in_123_type_cone": True,
        })

    owners, edges, _, neighbors = structure(k)
    equality_indices = [i for i, t in enumerate(types)
                        if all(abs(t[a]-t[b]) == 1 for a, b in edges)
                        and all(sum(t[u] == t[v]-1 for u in neighbors[v]) == 1
                                for v in range(k))]
    equality_types = [types[i] for i in equality_indices]

    # Independent direct description: the unique descending route exits through
    # one of the six leaves of the extended path.
    expected = {}
    for leaf, owner in enumerate(owners):
        core = tuple(abs(v-owner)-owner for v in range(k))
        leaves = tuple(core[v] + (-1 if j == leaf else 1)
                       for j, v in enumerate(owners))
        expected[core+leaves] = Fraction(1, 3*2**owner)
    assert set(equality_types) == set(expected)

    diagonal = [[power_two(-t[row]) for t in equality_types]
                for row in range(2*k+2)]
    _, pivots = exact_rref(diagonal)
    augmented, augmented_pivots = exact_rref([row+[Fraction(1)] for row in diagonal])
    assert len(pivots) == 6 and augmented_pivots == list(range(6))
    solution = [augmented[i][-1] for i in range(6)]
    assert solution == [expected[t] for t in equality_types]
    assert sum(solution) == 1
    assert all(sum(a*x for a, x in zip(row, solution)) == 1 for row in diagonal)

    return {
        "schema_version": 1,
        "status": "EXACT INTEGER/RATIONAL FINITE CHECKS PASSED; SEE ANALYTIC PROOF",
        "path_vertices": 4,
        "distance_types": len(types),
        "pair_projections": projection_results,
        "equality_support_size": len(equality_types),
        "diagonal_constraint_rank": len(pivots),
        "unique_normalized_masses": [
            {"type_index": index, "type": list(t), "mass": str(mass)}
            for index, t, mass in zip(equality_indices, equality_types, solution)],
        "scope": [
            "Checks exact local-type projections and equality-support algebra.",
            "Does not certify a numerical global maximum of the P4 cone.",
            "The strict-gap conclusion additionally uses the upstream pair theorem and the analytic arguments in GLOBAL_STRATEGY.md.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    rendered = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
