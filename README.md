# Hydrogen Orbitals

Computes hydrogen atom orbitals from the analytical solution of the time-independent Schrödinger equation
and saves the wave function on a spherical mesh as NumPy files. Two Python scripts then plot them:
2D slices (Matplotlib) and interactive 3D isosurfaces (Plotly).

```
orbitals.exe -n 4 -l 3 -m 1          →  results/n4_l3_m1/{r,theta,phi,psi,qn}.npy
python orbitals/visualize2d.py results/n4_l3_m1
python orbitals/visualize3d.py results/n4_l3_m1
```

## Physics and conventions

- **Units:** atomic units. Distances are in Bohr radii (a0 = 1). The proton is treated as infinitely heavy (no reduced-mass correction).
- **Wave function:** Ψn,l,m(r,θ,ϕ) = Rn,l(r) · Yl,m(θ,ϕ), with
  - Rn,l(r) = N · e^(−ρ) · (2ρ)^l · L^(2l+1)_(n−l−1)(2ρ), where ρ = r/n. L is the generalized Laguerre polynomial in the modern convention (L₁^(α)(x) = 1 + α − x). Rn,l is normalized: ∫ R² r² dr = 1.
  - Yl,m is the **complex** spherical harmonic, normalized to 1 on the sphere, **with the Condon–Shortley phase** (−1)^m. This is the convention of SciPy's `sph_harm_y` and C++17's `std::sph_legendre`.
- **Quantum numbers:** 1 ≤ n ≤ 20, 0 ≤ l < n, −l ≤ m ≤ l. The limit n ≤ 20 is `QuantumState::max_n`. It is the range over which the radial cut-off (`QuantumState::RadialExtent`, r_max = n(2n + 12)) was verified to leave less than 5·10⁻⁶ of the probability outside the mesh.
- **Mesh** (`Mesh::for_state`):
  - r ∈ [0, r_max] with a step of 0.2;
  - θ ∈ [0, π] with an odd number of nodes (≥ 21 and ≥ 2l + 1);
  - ϕ ∈ [0, 2π], with at least max(40, 4|m|) intervals; the last node repeats ϕ = 0.
- **Normalization check:** the program integrates |Ψ|² over the mesh and prints the total probability P. For θ it uses the Clenshaw–Curtis rule, which is exact for |Yl,m|²; for r and ϕ it uses the trapezoid rule. P ≈ 1 up to the radial quadrature error, which is about 1·10⁻⁴ for n = 1 and decreases with n.

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
orbitals.exe -n <n> -l <l> -m <m> [-o <dir>]
```

The output directory defaults to `results/n<n>_l<l>_m<m>`. It contains:

| file        | contents                                              |
|-------------|-------------------------------------------------------|
| `r.npy`     | radial nodes, shape (Nr,)                              |
| `theta.npy` | polar nodes θ, shape (Nθ,)                             |
| `phi.npy`   | azimuthal nodes ϕ, shape (Nϕ,)                         |
| `psi.npy`   | complex128 Ψ, shape (Nr, Nθ, Nϕ), C order              |
| `qn.npy`    | int32 quantum numbers [n, l, m]                        |

The size of `psi.npy` grows quickly with n: 5 MB for n = 4, and about 250 MB for n = 20, l = m = 19.

## Visualize

```sh
pip install -r requirements.txt

# 2D slice: --plane xz (through the Z axis, at azimuth --angle) or xy (cone θ = --angle, default 90°)
python orbitals/visualize2d.py results/n4_l3_m1 --quantity density|real|imag --plane xz --angle 0

# 3D isosurfaces enclosing the given fractions of ∫|Ψ|² (opens a browser, or writes --html)
python orbitals/visualize3d.py results/n4_l3_m1 --fraction 0.5 0.9 --grid 80
```

`--quantity real` and `imag` show Re Ψ and Im Ψ of the **complex** orbital. For m ≠ 0, these are the cos(mϕ) and sin(mϕ) combinations. They are not the chemists' real orbitals labelled by m (px, py, dxy …).

## Tests

```sh
ctest --test-dir build        # or: ctest --preset release
```

- `test_orbitals` (C++) checks:
  - the Laguerre and Legendre functions, against closed forms and against `std::assoc_laguerre` / `std::sph_legendre` where the standard library provides them (MSVC STL, libstdc++);
  - that the spherical harmonics are orthonormal, for both the complex and the real basis;
  - that the tail beyond `RadialExtent()` is < 5·10⁻⁶ for every (n, l);
  - that the mesh probability is ≈ 1.
- `scipy_*` tests run the program for a few states and compare `psi.npy` with an independent SciPy formula (`tests/check_against_scipy.py`). They are added only if Python with NumPy and SciPy is found when CMake configures.

## Layout

```
orbitals/main.cpp             command line (CLI11), runs the computation, saves the mesh
orbitals/quantum_state.hpp    QuantumState: Rn,l, Yl,m, Ψ, radial extent
orbitals/mesh.hpp             Mesh: the grid, Ψ on the grid, quadrature, .npy output
orbitals/math/                Laguerre and Legendre polynomials, spherical harmonics, factorial, Kahan sum
orbitals/npy.hpp              minimal .npy writer
orbitals/visualize2d.py       2D slices
orbitals/visualize3d.py       3D isosurfaces
tests/                        C++ unit tests and the SciPy cross-check
```

The original program was written by Derek Harrison.
