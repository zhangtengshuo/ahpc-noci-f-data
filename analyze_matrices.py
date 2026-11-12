#!/usr/bin/env python3
"""Calculate energies and pairwise couplings from a distributed NOCI result."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

EH_TO_EV = 27.211386245988


def read_result(path):
    if path.suffix == ".json":
        result = json.loads(path.read_text())
        return (result["labels"], np.asarray(result["normalized"]["hamiltonian"], dtype=float),
                np.asarray(result["normalized"]["overlap"], dtype=float),
                result.get("baseline_caspt2", {}).get("shifted_hamiltonian"))
    matrices = {}
    labels = None
    for line in path.read_text().splitlines():
        fields = line.split("|")
        if fields[0] == "labels":
            labels = fields[1].split(",")
        elif fields[0] in ("norm_hamiltonian", "norm_overlap"):
            matrices.setdefault(fields[0], []).append([float(x) for x in fields[1].split()])
    return labels, np.asarray(matrices["norm_hamiltonian"]), np.asarray(matrices["norm_overlap"]), None


def write(path, columns, rows):
    with path.open("w", newline="") as stream:
        out = csv.writer(stream)
        out.writerow(columns)
        out.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--flip-state", action="append", default=[], help="State phase to flip in the exported matrices/couplings.")
    args = parser.parse_args()
    labels, h, s, shifted = read_result(args.result)
    for label in args.flip_state:
        if label not in labels:
            parser.error(f"Unknown state: {label}")
    phases = np.array([-1.0 if label in args.flip_state else 1.0 for label in labels])
    transform = phases[:, None] * phases[None, :]
    s = s * transform
    eigenvalues, vectors = np.linalg.eigh(s)
    if np.min(eigenvalues) <= 0:
        raise ValueError("The overlap matrix must be positive definite.")
    inverse_root = (vectors * eigenvalues**-0.5) @ vectors.T
    energies, couplings = [], []
    versions = [("unshifted", h)]
    if shifted is not None:
        versions.append(("CASPT2_shifted", np.asarray(shifted)))
    for name, matrix in versions:
        matrix = matrix * transform
        roots = np.linalg.eigvalsh(inverse_root @ matrix @ inverse_root)
        for i, label in enumerate(labels):
            energies.append((name, "diabatic", label, matrix[i, i], (matrix[i, i] - matrix[0, 0]) * EH_TO_EV))
        for i, energy in enumerate(roots, 1):
            energies.append((name, "adiabatic", str(i), energy, (energy - roots[0]) * EH_TO_EV))
        for i in range(len(labels)):
            for j in range(i):
                coupling = (matrix[i, j] - (matrix[i, i] + matrix[j, j]) * s[i, j] / 2) / (1 - s[i, j]**2)
                couplings.append((name, labels[i], labels[j], matrix[i, j], s[i, j], coupling, coupling * EH_TO_EV * 1000))
    args.output_directory.mkdir(parents=True, exist_ok=True)
    write(args.output_directory / "energies.csv", ["hamiltonian", "kind", "state_or_root", "energy_Eh", "relative_energy_eV"], energies)
    write(args.output_directory / "couplings.csv", ["hamiltonian", "state_i", "state_j", "Hij_Eh", "Sij", "coupling_Eh", "coupling_meV"], couplings)


if __name__ == "__main__":
    main()
