# Hydrogen Orbitals

Computes hydrogen atom orbitals from the analytical solution of the time-independent Schrödinger equation
and saves the wave function on a spherical mesh as NumPy files. Two Python scripts then plot them:
2D slices (Matplotlib) and interactive 3D isosurfaces (Plotly).

This is a fork of Derek Harrison's [hydrogen-orbitals](https://github.com/derekharrison/hydrogen-orbitals); see [Origin](#origin).

```
orbitals.exe -n 4 -l 3 -m 1               →  results/n4_l3_m1/{r,theta,phi,psi,qn}.npy, info.json
orbitals.exe -n 3 -l 2 -m -2 --basis real  →  results/n3_l2_m-2_real/   (the 3d_xy orbital)
python orbitals/visualize2d.py results/n4_l3_m1
python orbitals/visualize3d.py results/n3_l2_m-2_real --quantity real
```

## Physics and conventions

- **Units:** atomic units.

  Distances are in Bohr radii ($a_0 = 1$). The proton is treated as infinitely heavy (no reduced-mass correction).

- **Wave function:** $Ψ_{n,l,m}(r,θ,ϕ) = R_{n,l}(r) · Y_{l,m}(θ,ϕ)$, with

  - $R_{n,l}(r) = N · e^{−ρ} · (2ρ)^l · L^{(2l+1)}_{n−l−1}(2ρ)$  is the radial part of $Ψ$, where:
  
    - $ρ = r/(n·a_0)$,
    - $N = \sqrt{ \left( \frac{2}{n·a_0} \right)^3 · \frac{(n-l-1)!}{2n·(n+l)!} }$ - normalizing factor.
    - $L$ is the generalized Laguerre polynomial in the modern convention: $L_1^{(α)}(x) = 1 + α − x$.
    
    $R_{n,l}$ is normalized: $\int_0^\infty R_{n,l}^2 \, r^2 \, dr = 1$.

  - $Y_{l,m}$ for `--basis complex` (which is the default) is the **complex** spherical harmonic, normalized to 1 on the sphere, **with the Condon–Shortley phase** $(−1)^m$ (this is the convention of SciPy's `sph_harm_y` and C++17's `std::sph_legendre`). These states are eigenstates of $L_z$ with eigenvalue $m$.

  - $Y_{l,m}$ for `--basis real` is the **real** spherical harmonic:

    - $\sqrt{2}·(−1)^m·\mathrm{Re}(Y_{l,m})$ for $m > 0$,
    - $Y_{l,0}$ for $m = 0$,
    - $\sqrt{2}·(−1)^m·\mathrm{Im}(Y_{l,|m|})$ for $m < 0$.
    
    These are the orbitals of chemistry, which are not eigenstates of $L_z$ for $m ≠ 0$. Each is proportional to a Cartesian polynomial divided by $r^l$, and is named after it:
    
    - for $l = 1, m = 1, −1, 0$ give $p_x$, $p_y$, $p_z$;
    - for $l = 2, m = −2 … 2$ give $d_{xy}$, $d_{yz}$, $d_{z^2}$, $d_{xz}$, $d_{x^2-y^2}$;
    - for $l = 3$ gives the seven $f$ orbitals.
    
    The positive lobe points along the positive axes, so $p_x > 0$ for $x > 0$.

- **Quantum numbers:** $1 ≤ n ≤ 20$, $0 ≤ l < n$, $−l ≤ m ≤ l$.

  The limit $n ≤ 20$ (`QuantumState::max_n`) is the range over which the radial cut-off (`QuantumState::RadialExtent`, $r_{max} = n·(2n + 12)$) was verified to leave less than $5·10^{-6}$ of the probability outside the mesh.

- **Mesh** (`Mesh::for_state`; the defaults can be changed with the options under [Run](#run)):

  - $r ∈ [0, r_{max}]$ with a step of 0.2;

  - $θ ∈ [0, π]$ with an odd number of nodes (default $max(21, 2l + 1)$).
  
    The node count must be odd, so that there is a node at $θ = π/2$, and $≥ 2l + 1$, so that the $θ$ integral is exact.

  - $ϕ ∈ [0, 2π]$, with $max(40, 4|m|)$ intervals by default; the last node repeats $ϕ = 0$.

    The number of intervals must be even, so that every half-plane has its opposite in the $xz$ plots, and $> 2|m|$, so that the $ϕ$ integral is exact for real orbitals.

- **Checks:** the program integrates over the mesh and prints a table comparing each result with the exact formula. It also writes the same numbers to `info.json`:

  | quantity              | exact value                                                       |
  |-----------------------|-------------------------------------------------------------------|
  | total probability $P$ | $1$                                                               |
  | $⟨r⟩$                 | $(3n² − l(l+1)) / 2$                                              |
  | $⟨r²⟩$                | $n² (5n² + 1 − 3l(l+1)) / 2$                                      |
  | $⟨1/r⟩$               | $1 / n²$                                                          |
  | energy $E = −⟨1/r⟩/2$ | $−1 / (2n²)$ $\mathrm{E_h}$ (1 $\mathrm{E_h}$ = 27.211386 eV)     |

  The energy on the mesh comes from the virial theorem ($⟨V⟩ = 2E$ for the Coulomb potential), so it tests $⟨1/r⟩$ and the radial quadrature, not a separate solver.
  - **Quadrature:**

    - for $θ$, the *Clenshaw–Curtis* rule, which is exact for $|Y_{l,m}|²$.
    - for $ϕ$, the periodic trapezoid rule, which is exact here.
    - for $r$, the *Gregory* rule (`Mesh::gregory_weights`): the trapezoid rule with 6 correction terms at each end of the range. Its weights are all positive, and it is exact for polynomials up to degree 6.

  - **Why not the plain trapezoid rule:** for $s$ states ($l = 0$), the $⟨1/r⟩$ integrand $r·R²$ rises linearly from $r = 0$, which is the trapezoid rule's worst case. It made $⟨1/r⟩$ and $E$ 1.3% off for $1s$; the Gregory rule brings that to $4·10^{-5}$, and it converges much faster as $h$ shrinks (about $h^8$ instead of $h^2$).

  - **Typical accuracy, default step 0.2:** every check is within about $1·10^{-4}$ for $n = 1$ (the worst case: only 5 steps per decay length of the $1s$ wave function), and better for larger $n$. For large $n$ with $l=0$, the remaining deficit in $P$ (about $5·10^{-6}$ for $n = 20$) is the probability beyond $r_{max}$, not integration error. A smaller `--radial-step` improves $n = 1$ quickly: 0.1 gives better than $1·10^{-6}$.

## Build

Requirements: CMake ≥ 3.24, a C++17 or later compiler (C++20 by default, `-DCMAKE_CXX_STANDARD=17` for C++17; tested with clang++ 23 on Windows), Ninja.
[CLI11](https://github.com/CLIUtils/CLI11) ≥ 2.7 is found with `find_package`, or downloaded automatically.

```sh
cmake --preset release        # configures ./build (Ninja, Release); see CMakePresets.json
cmake --build build
```

To set the compiler or a local CLI11 path for your machine only, put a `CMakeUserPresets.json` next to `CMakePresets.json` (it is git-ignored), for example:

```json
{
  "version": 3,
  "configurePresets": [{
    "name": "local-release", "inherits": "release",
    "cacheVariables": { "CMAKE_CXX_COMPILER": "clang++", "CMAKE_PREFIX_PATH": "path/to/CLI11" }
  }]
}
```

## Run

```
orbitals.exe -n <n> -l <l> -m <m> [--basis complex|real] [--radial-step <h>]
             [--polar-points <Nθ>] [--azimuthal-points <Nϕ>] [-o <dir>]
```

| option               | default                                 | constraint                                                                |
|----------------------|-----------------------------------------|---------------------------------------------------------------------------|
| `--basis`            | `complex`                               | `complex` ($e^{imϕ}$, eigenstates of $L_z$) or `real` ($p_x$, $d_{xy}$ …) |
| `--radial-step`      | $0.2$ (Bohr radii)                      | $0 < h ≤ r_{max} / 2$                                                    |
| `--polar-points`     | $max(21, 2l + 1)$                      | odd, $≥ max(3, 2l + 1)$                                                  |
| `--azimuthal-points` | $max(40, 4 \lvert m \rvert) + 1$       | odd, $≥ 5$, $> 2 \lvert m \rvert + 1$ (the node at $2π$ included)         |

A mesh whose `psi` would exceed 4 GB is refused. Invalid options are reported, and the program exits with a non-zero code.

The output directory defaults to `results/n<n>_l<l>_m<m>`, with `_real` appended for the real basis. It contains:

| file        | contents                                              |
|-------------|-------------------------------------------------------|
| `r.npy`     | radial nodes, shape (Nr,)                              |
| `theta.npy` | polar nodes θ, shape (Nθ,)                             |
| `phi.npy`   | azimuthal nodes ϕ, shape (Nϕ,)                         |
| `psi.npy`   | complex128 Ψ, shape (Nr, Nθ, Nϕ), C order              |
| `qn.npy`    | int32 quantum numbers [n, l, m]                        |
| `info.json` | basis, orbital name, mesh shape and step, and the checks (each as `{"grid": …, "exact": …}`) |

For the real basis, `psi.npy` is still complex128, with a zero imaginary part, so all the tools read both bases the same way.
The size of `psi.npy` grows quickly with $n$: 5 MB for $n = 4$, and about 250 MB for $n = 20$, $l = m = 19$. Finer angular meshes multiply that: `--polar-points 91 --azimuthal-points 181` is about 20 times the default.

## Visualize

```sh
pip install -r requirements.txt

# 2D slice: --plane xz (through the Z axis, at azimuth --angle) or xy (cone θ = --angle, default 90°)
python orbitals/visualize2d.py results/n4_l3_m1 --quantity density|real|imag --plane xz --angle 0

# 3D isosurfaces enclosing the given fractions of ∫|Ψ|² (opens a browser, or writes --html)
python orbitals/visualize3d.py results/n4_l3_m1 --fraction 0.5 0.9 --grid 80
```

`--quantity real` and `imag` show $\mathrm{Re}(Ψ)$ and $\mathrm{Im}(Ψ)$. For a complex orbital with $m ≠ 0$, these are the $\cos(mϕ)$ and $\sin(mϕ)$ combinations. To show a chemists' real orbital ($p_y$, $d_{xy}$ …) with its signed lobes, compute it with `--basis real` and plot `--quantity real`; its imaginary part is zero. The plot titles give the orbital name from `info.json`.

## Tests

```sh
ctest --test-dir build        # or: ctest --preset release
```

- `test_orbitals` (C++) checks:
  - the Laguerre and Legendre functions, against closed forms and against `std::assoc_laguerre` / `std::sph_legendre` where the standard library provides them (MSVC STL, libstdc++);
  - that the spherical harmonics are orthonormal, for both the complex and the real basis;
  - that the tail beyond `RadialExtent()` is < 5·10⁻⁶ for every (n, l);
  - that the mesh probability, ⟨r⟩, ⟨r²⟩ and ⟨1/r⟩ agree with the exact values, in both bases;
  - that the Gregory weights are exact for polynomials up to degree 6, and that the radial integration converges at high order (the error falls about 200× per halving of the step);
  - that invalid resolutions are rejected;
  - the orbital names.
- `cli_*` tests check that invalid command-line options make the program fail.
- `scipy_*` tests run the program for a few states in both bases and compare `psi.npy` with an independent SciPy formula (`tests/check_against_scipy.py`). They are added only if Python with NumPy and SciPy is found when CMake configures.

## Layout

```
orbitals/main.cpp             command line (CLI11), runs the computation, saves the mesh
orbitals/quantum_state.hpp    QuantumState: Rn,l, Yl,m (complex or real basis), Ψ, radial extent, orbital name, exact values
orbitals/mesh.hpp             Mesh: the grid and its resolution, Ψ on the grid, quadrature and checks, .npy and info.json output
orbitals/math/                Laguerre and Legendre polynomials, spherical harmonics, factorial, Kahan sum
orbitals/npy.hpp              minimal .npy writer
orbitals/visualize2d.py       2D slices
orbitals/visualize3d.py       3D isosurfaces
tests/                        C++ unit tests and the SciPy cross-check
```

## Origin

This repository is a fork of [derekharrison/hydrogen-orbitals](https://github.com/derekharrison/hydrogen-orbitals), written by Derek Harrison in 2021. It contained a C++ program that evaluates the orbitals from the analytical solution on a spherical grid, and a Python script that plots them.

Its history is kept here: the first 36 commits are the original repository's. In 2024 the original files were moved aside and reworked step by step. They have since been replaced by the code described above, which keeps the original approach but is largely rewritten:

- header-only C++20 (or C++17) code, built with CMake;
- a command line (CLI11);
- NumPy output, with 2D and 3D visualizers;
- the real orbitals;
- numerical checks against exact values;
- tests.

The original repository does not specify a license, and this fork does not add one.
