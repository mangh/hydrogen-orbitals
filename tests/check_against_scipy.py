# ASCII output only: the Windows console may not be able to print Greek letters.
"""Compares psi saved by orbitals.exe with an independent formula built from SciPy.

    psi(r, theta, phi) = R_nl(r) * Y_lm(theta, phi)   (atomic units, a0 = 1)

    R_nl(r) = sqrt((2/n)^3 (n-l-1)! / (2n (n+l)!)) * exp(-r/n) * (2r/n)^l * L_(n-l-1)^(2l+1)(2r/n)
    Y_lm    = scipy.special.sph_harm_y (SciPy >= 1.15) or sph_harm (older), both with the Condon-Shortley phase

Exits with status 1 when the largest deviation, relative to max |psi|, exceeds the tolerance.

Usage:
    python check_against_scipy.py results/n4_l3_m1 [--tol 1e-10]
"""

import argparse
import sys
from math import factorial
from pathlib import Path

import numpy as np
from scipy import special


def spherical_harmonic(l, m, theta, phi):
    if hasattr(special, "sph_harm_y"):
        return special.sph_harm_y(l, m, theta, phi)
    return special.sph_harm(m, l, phi, theta)  # older SciPy: (order, degree, azimuth, polar)


def reference_psi(n, l, m, r, theta, phi):
    rho = 2.0 * r / n
    norm = np.sqrt((2.0 / n) ** 3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
    radial = norm * np.exp(-rho / 2) * rho ** l * special.eval_genlaguerre(n - l - 1, 2 * l + 1, rho)
    T, P = np.meshgrid(theta, phi, indexing="ij")
    angular = spherical_harmonic(l, m, T, P)
    return radial[:, None, None] * angular[None, :, :]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data", type=Path, help="directory written by orbitals.exe")
    parser.add_argument("--tol", type=float, default=1e-10, help="largest allowed |psi - reference| / max|reference|")
    args = parser.parse_args()

    r, theta, phi, psi = (np.load(args.data / f"{name}.npy") for name in ("r", "theta", "phi", "psi"))
    n, l, m = (int(q) for q in np.load(args.data / "qn.npy"))

    expected_shape = (len(r), len(theta), len(phi))
    if psi.shape != expected_shape:
        sys.exit(f"psi has shape {psi.shape}, expected {expected_shape}")

    ref = reference_psi(n, l, m, r, theta, phi)
    error = np.abs(psi - ref).max() / np.abs(ref).max()
    ok = error <= args.tol
    print(f"n={n}, l={l}, m={m}: max |psi - reference| / max |reference| = {error:.3g}"
          f" ({'ok' if ok else 'FAILED'}, tolerance {args.tol:g})")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
