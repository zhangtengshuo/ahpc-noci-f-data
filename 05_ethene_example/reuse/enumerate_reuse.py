#!/usr/bin/env python3
"""Enumerate the complete occupation-pattern example of SI Fig. S1."""
import argparse
import csv
from collections import Counter
from itertools import product
from pathlib import Path


def write(path, columns, rows):
    with path.open("w", newline="") as stream:
        out = csv.writer(stream)
        out.writerow(columns)
        out.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    output = args.output_directory
    output.mkdir(parents=True, exist_ok=True)
    # Each neutral MS=0 monomer has one alpha and one beta electron in two orbitals.
    left_strings = ["1010", "1001", "0110", "0101"]
    left = list(product(left_strings, repeat=2))
    # TT coupled to a singlet contains MS=(+1,-1), (-1,+1), and four (0,0) terms.
    right_alpha = ["1100", "0011", "1010", "1001", "0110", "0101"]
    right = [(a, "".join("1" if x == "0" else "0" for x in a)) for a in right_alpha]
    identifiers = {(a, b): i for i, (a, b) in enumerate(product(left_strings, right_alpha), 1)}
    alpha, beta = Counter(), Counter()
    mappings = []
    for i, (la, lb) in enumerate(left, 1):
        for j, (ra, rb) in enumerate(right, 1):
            alpha[(la, ra)] += 1
            beta[(lb, rb)] += 1
            mappings.append((i, j, identifiers[(la, ra)], identifiers[(lb, rb)]))
    assert len(mappings) == 96 and sum(alpha.values()) + sum(beta.values()) == 192
    assert set(alpha) == set(beta) and len(alpha.keys() | beta.keys()) == 24
    determinants = [(state, i, a, b) for state, dets in [("S1S0", left), ("T1T1", right)]
                    for i, (a, b) in enumerate(dets, 1)]
    write(output / "determinants.csv", ["state", "determinant", "alpha", "beta"], determinants)
    write(output / "pair_mapping.csv", ["left_determinant", "right_determinant", "alpha_pair", "beta_pair"], mappings)
    write(output / "unique_pairs.csv", ["pair", "left_string", "right_string", "alpha_occurrences", "beta_occurrences"],
          [(index, a, b, alpha[(a, b)], beta[(a, b)]) for (a, b), index in identifiers.items()])
    print("16 x 6 = 96 determinant pairs; 192 single-spin occurrences; 24 unique pairs.")


if __name__ == "__main__":
    main()
