# The docstring (also the --help text) is kept ASCII, so that it prints on any Windows console code page.
"""Shows a hydrogen orbital saved by orbitals.exe (see Mesh::save) as interactive 3D isosurfaces.

The wave function psi[r, theta, phi] is resampled from the spherical mesh onto a Cartesian cube
and drawn with Plotly in the browser (rotate / zoom with the mouse). Each surface is chosen so that
it encloses a given fraction of the integral of |psi|^2 (or of (Re psi)^2, (Im psi)^2 for the real
and imaginary parts, whose positive and negative lobes are drawn in red and blue).

Requires NumPy, SciPy and Plotly: pip install -r requirements.txt

Examples:
    python visualize3d.py results/n4_l3_m1
    python visualize3d.py results/n3_l2_m0 --quantity real
    python visualize3d.py results/n4_l3_m1 --fraction 0.5 0.9 --grid 80 --html n4_l3_m1.html
"""

import argparse
import sys
from pathlib import Path

try:
    import numpy as np
    import plotly.graph_objects as go
    from plotly.colors import sample_colorscale
    from scipy.interpolate import RegularGridInterpolator
except ImportError as e:
    sys.exit(f"visualize3d.py needs {e.name}: pip install -r requirements.txt")

from visualize2d import QUANTITIES

POSITIVE, NEGATIVE = "#d6604d", "#4393c3"  # lobe colors for Re Ψ / Im Ψ (as in RdBu)


def radial_extent(r, theta, phi, psi, fraction):
    """Radius of the ball that holds the given fraction of the probability ∫|Ψ|² dV."""
    # ϕ: the last node duplicates the first, so it is left out of the sum.
    density = (np.abs(psi[:, :, :-1]) ** 2).sum(axis=2)  # (r, θ)
    radial = (density * np.sin(theta)).sum(axis=1) * r ** 2  # ∝ P(r)
    cumulative = np.cumsum(radial)
    return r[min(np.searchsorted(cumulative, fraction * cumulative[-1]), len(r) - 1)]


def resample(r, theta, phi, psi, half_width, points, method):
    """Ψ on a Cartesian grid of points³ nodes spanning [-half_width, half_width]³."""
    axis = np.linspace(-half_width, half_width, points)
    X, Y, Z = np.meshgrid(axis, axis, axis, indexing="ij")
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    T = np.arctan2(np.hypot(X, Y), Z)  # θ ∈ [0, π], well defined at the origin
    P = np.mod(np.arctan2(Y, X), 2 * np.pi)  # ϕ ∈ [0, 2π), covered by the ϕ grid (it ends at 2π)

    values = np.stack([psi.real, psi.imag], axis=-1)
    interp = RegularGridInterpolator((r, theta, phi), values, method=method, bounds_error=False, fill_value=0.0)
    re_im = interp(np.stack([R, T, P], axis=-1))  # beyond r_max (cube corners): 0
    return X, Y, Z, re_im[..., 0] + 1j * re_im[..., 1]


def enclosing_level(weight, fraction):
    """Threshold c such that the nodes with weight > c hold the given fraction of the total weight."""
    w = np.sort(weight.ravel())[::-1]
    cumulative = np.cumsum(w)
    return w[min(np.searchsorted(cumulative, fraction * cumulative[-1]), len(w) - 1)]


def surface(X, Y, Z, F, level, opacity, color=None):
    """Isosurface F = level in a uniform color, or (color=None) the pair F = ±level in red / blue."""
    lo, hi, count, colors = (level, level, 1, [color, color]) if color else (-level, level, 2, [NEGATIVE, POSITIVE])
    # float32 halves the size of the data embedded in the figure (and in the HTML file).
    return go.Isosurface(
        x=X.ravel().astype(np.float32), y=Y.ravel().astype(np.float32), z=Z.ravel().astype(np.float32),
        value=F.ravel().astype(np.float32),
        isomin=lo, isomax=hi, surface_count=count, cmin=lo, cmax=hi,
        colorscale=[[0, colors[0]], [1, colors[1]]], showscale=False,
        caps=dict(x_show=False, y_show=False, z_show=False),
        opacity=opacity, hoverinfo="skip",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data", type=Path, help="directory written by orbitals.exe")
    parser.add_argument("--quantity", choices=QUANTITIES, default="density")
    parser.add_argument("--fraction", type=float, nargs="+", default=[0.9],
                        help="share(s) of the integral enclosed by the surface(s), e.g. 0.5 0.9 for nested shells")
    parser.add_argument("--grid", type=int, default=64, help="Cartesian nodes per axis")
    parser.add_argument("--extent", type=float, default=0.995,
                        help="the cube encloses the ball holding this share of the probability")
    parser.add_argument("--method", choices=["linear", "cubic"], default="cubic",
                        help="interpolation from the spherical mesh to the Cartesian grid (linear: faster, more faceted)")
    parser.add_argument("--html", type=Path, help="write the figure to this HTML file instead of opening the browser")
    args = parser.parse_args()

    if not all(0 < p < 1 for p in args.fraction):
        parser.error("--fraction values must lie between 0 and 1")

    r, theta, phi, psi = (np.load(args.data / f"{name}.npy") for name in ("r", "theta", "phi", "psi"))
    n, l, m = np.load(args.data / "qn.npy")

    half_width = radial_extent(r, theta, phi, psi, args.extent)
    X, Y, Z, psi_xyz = resample(r, theta, phi, psi, half_width, args.grid, args.method)

    func, label = QUANTITIES[args.quantity]
    F = func(psi_xyz)
    signed = args.quantity != "density"
    weight = F ** 2 if signed else F
    if not weight.any():
        sys.exit(f"{args.quantity} is zero everywhere (n={n}, l={l}, m={m}), nothing to draw")

    traces = []
    # Inner surfaces (small fractions) opaque, the outer ones translucent so that the inner stay visible.
    fractions = sorted(set(args.fraction))
    for i, fraction in enumerate(fractions):
        level = enclosing_level(weight, fraction)
        opacity = 1.0 if i == 0 else 0.3
        if signed:
            level = np.sqrt(level)
            traces.append(surface(X, Y, Z, F, level, opacity))
        else:
            shade = 0.8 - 0.5 * i / max(len(fractions) - 1, 1)  # inner bright, outer dark
            traces.append(surface(X, Y, Z, F, level, opacity, sample_colorscale("Inferno", [shade])[0]))
        # ASCII only: the Windows console may not be able to print Ψ.
        print(f"{fraction:.0%} surface: {args.quantity} = {'+/-' if signed else ''}{level:.4g}")

    enclosed = ", ".join(f"{p:.0%}" for p in fractions)
    integrand = "|Ψ|²" if not signed else f"({label})²"
    fig = go.Figure(traces)
    axis = dict(range=[-half_width, half_width])
    fig.update_layout(
        title=f"{label}   n={n}, l={l}, m={m}   (surfaces enclosing {enclosed} of ∫{integrand})",
        scene=dict(xaxis=dict(title="X", **axis), yaxis=dict(title="Y", **axis), zaxis=dict(title="Z", **axis),
                   aspectmode="cube"),
        margin=dict(l=0, r=0, b=0, t=40),
    )

    if args.html:
        fig.write_html(args.html, include_plotlyjs="cdn")
        print(f"Saved {args.html.resolve()}")
    else:
        fig.show()


if __name__ == "__main__":
    main()
