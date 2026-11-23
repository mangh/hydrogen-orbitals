"""Plots a planar slice of a hydrogen orbital saved by orbitals.exe (see Mesh::save).

The data directory holds the mesh coordinates r, theta, phi and the complex wave function
psi[r, theta, phi] as NumPy arrays; any plane and quantity can be cut out of it here.

Examples:
    python visualize2d.py results/n4_l3_m1
    python visualize2d.py results/n4_l3_m1 --quantity real --plane xy
    python visualize2d.py results/n4_l3_m1 --plane xz --angle 45
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

QUANTITIES = {
    "density": (lambda psi: np.abs(psi) ** 2, "|Ψ|²"),
    "real": (lambda psi: psi.real, "Re Ψ"),
    "imag": (lambda psi: psi.imag, "Im Ψ"),
}


def nearest(grid, angle):
    """Index of the grid node closest to the given angle (mod 2π)."""
    diff = np.angle(np.exp(1j * (grid - angle)))
    return int(np.argmin(np.abs(diff)))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data", type=Path, help="directory written by orbitals.exe")
    parser.add_argument("--quantity", choices=QUANTITIES, default="density")
    parser.add_argument("--plane", choices=["xz", "xy"], default="xz",
                        help="xz: vertical plane through the Z axis; xy: plane (or cone) of constant theta")
    parser.add_argument("--angle", type=float,
                        help="xz: azimuth phi of the plane in degrees (default 0); xy: polar angle theta in degrees (default 90)")
    parser.add_argument("--levels", type=int, default=100, help="number of contour levels")
    args = parser.parse_args()

    r, theta, phi, psi = (np.load(args.data / f"{name}.npy") for name in ("r", "theta", "phi", "psi"))
    n, l, m = np.load(args.data / "qn.npy")

    func, label = QUANTITIES[args.quantity]
    f = func(psi)  # shape (r, θ, ϕ)

    if args.plane == "xz":
        # Two opposite half-planes ϕ0 and ϕ0 + π, joined along θ into one full turn
        # (angle measured from +Z: 0..π on the ϕ0 side, π..2π on the opposite side).
        phi0 = np.radians(args.angle if args.angle is not None else 0.0)
        i0, i1 = nearest(phi, phi0), nearest(phi, phi0 + np.pi)
        angle = np.concatenate([theta, 2 * np.pi - theta[::-1]])
        F = np.concatenate([f[:, :, i0], f[:, ::-1, i1]], axis=1)
        R, A = np.meshgrid(r, angle, indexing="ij")
        X, Y = R * np.sin(A), R * np.cos(A)
        xlabel, ylabel = ("X" if i0 == 0 else f"horizontal axis at ϕ = {np.degrees(phi[i0]):.0f}°"), "Z"
        where = f"ϕ = {np.degrees(phi[i0]):.0f}°"
    else:
        theta0 = np.radians(args.angle if args.angle is not None else 90.0)
        j = int(np.argmin(np.abs(theta - theta0)))
        if np.sin(theta[j]) < 1e-9:
            parser.error(f"--angle {args.angle}: the cone theta = {np.degrees(theta[j]):.0f} deg degenerates to the Z axis")  # ASCII: the console may not print Greek letters
        F = f[:, j, :]
        # For θ ≠ 90° the nodes lie on a cone; it is drawn as its projection onto the XY-plane (radius r·sinθ).
        R, P = np.meshgrid(r * np.sin(theta[j]), phi, indexing="ij")
        X, Y = R * np.cos(P), R * np.sin(P)
        flat = np.isclose(theta[j], np.pi / 2)
        xlabel, ylabel = ("X", "Y") if flat else ("X (projection of the cone)", "Y (projection of the cone)")
        where = f"θ = {np.degrees(theta[j]):.0f}°"

    vmax = np.abs(F).max() or 1.0
    if args.quantity == "density":
        levels, cmap = np.linspace(0, vmax, args.levels), "inferno"
    else:
        levels, cmap = np.linspace(-vmax, vmax, args.levels), "RdBu_r"

    fig, ax = plt.subplots()
    cont = ax.contourf(X, Y, F, levels=levels, cmap=cmap)
    fig.colorbar(cont, ax=ax, label=label)
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{label}   n={n}, l={l}, m={m}   ({args.plane.upper()}, {where})")
    plt.show()


if __name__ == "__main__":
    main()
