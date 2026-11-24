/*
 * test_orbitals.cpp
 *
 * Unit tests of the special functions, the wave function normalization and the mesh quadrature.
 *
 *   test_orbitals            - runs all tests,
 *   test_orbitals <name>...  - runs the named tests (this is how CTest calls it, see tests/CMakeLists.txt).
 *
 * The exit code is EXIT_FAILURE when any check fails.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "math/assoc_laguerre.hpp"
#include "math/assoc_legendre.hpp"
#include "math/constants.hpp"
#include "math/spherical_harmonic.hpp"
#include "mesh.hpp"
#include "numbers.hpp"
#include "quantum_state.hpp"

using namespace mangh::math;
using namespace mangh::hydrogenium;

namespace
{
    constexpr real_t PI = PI_V<real_t>;

    int g_failures = 0;

    /// Reports a failed check (printf-style message), and counts it.
    template<typename... Args>
    void check(bool ok, const char* format, Args... args)
    {
        if (ok) return;
        ++g_failures;
        std::printf("  FAILED: ");
        if constexpr (sizeof...(Args) == 0)
            std::fputs(format, stdout);
        else
            std::printf(format, args...);
        std::printf("\n");
    }

    /// |a - b| <= tol · max(1, |b|): relative for large values, absolute for small ones.
    bool close(real_t a, real_t b, real_t tol)
    {
        return std::abs(a - b) <= tol * std::max<real_t>(1.0, std::abs(b));
    }

    ///////////////////////////////////////////////////////////////////
    ///
    ///     assoc_laguerre
    ///

    void test_laguerre()
    {
        // Closed forms of the first three polynomials, for every α used by Rn,l (α = 2l+1 <= 39)
        for (unsigned alpha = 0; alpha < 2 * QuantumState::max_n; ++alpha)
        {
            const real_t a = alpha;
            for (real_t x : { 0.0, 0.3, 1.7, 5.0, 12.5 })
            {
                check(close(assoc_laguerre(alpha, 0, x), 1.0, 1e-14), "L(a=%u, n=0, x=%g)", alpha, x);
                check(close(assoc_laguerre(alpha, 1, x), 1.0 + a - x, 1e-14), "L(a=%u, n=1, x=%g)", alpha, x);
                const real_t L2 = x * x / 2.0 - (a + 2.0) * x + (a + 2.0) * (a + 1.0) / 2.0;
                check(close(assoc_laguerre(alpha, 2, x), L2, 1e-13), "L(a=%u, n=2, x=%g)", alpha, x);
            }
        }

#if defined(__cpp_lib_math_special_functions)
        // The C++17 special functions (MSVC STL, libstdc++; not libc++): std::assoc_laguerre(n, m, x) = L_n^(m)(x)
        for (unsigned alpha = 0; alpha < 2 * QuantumState::max_n; ++alpha)
            for (unsigned n = 0; n < QuantumState::max_n; ++n)
                for (real_t x : { 0.1, 1.0, 4.0, 10.0, 40.0 })
                {
                    const real_t expected = std::assoc_laguerre(n, alpha, x);
                    check(close(assoc_laguerre(alpha, n, x), expected, 1e-9),
                          "L(a=%u, n=%u, x=%g) = %.17g, std::assoc_laguerre = %.17g", alpha, n, x,
                          assoc_laguerre(alpha, n, x), expected);
                }
#else
        std::printf("  (std::assoc_laguerre not available: compared with the closed forms only)\n");
#endif
    }

    ///////////////////////////////////////////////////////////////////
    ///
    ///     ylm_assoc_legendre
    ///

    void test_legendre()
    {
        // Closed forms (with the Condon-Shortley phase): P~l,m(cosθ) = Yl,m(θ, 0)
        for (real_t theta : { 0.0, 0.4, 1.0, PI / 2, 2.5, PI })
        {
            const real_t c = std::cos(theta), s = std::sin(theta);
            check(close(ylm_assoc_legendre(0u, 0u, c), std::sqrt(1.0 / (4 * PI)), 1e-15), "Y00(%g)", theta);
            check(close(ylm_assoc_legendre(1u, 0u, c), std::sqrt(3.0 / (4 * PI)) * c, 1e-15), "Y10(%g)", theta);
            check(close(ylm_assoc_legendre(1u, 1u, c), -std::sqrt(3.0 / (8 * PI)) * s, 1e-15), "Y11(%g)", theta);
            check(close(ylm_assoc_legendre(2u, 0u, c), std::sqrt(5.0 / (16 * PI)) * (3 * c * c - 1), 1e-15), "Y20(%g)", theta);
            check(close(ylm_assoc_legendre(2u, 1u, c), -std::sqrt(15.0 / (8 * PI)) * s * c, 1e-15), "Y21(%g)", theta);
            check(close(ylm_assoc_legendre(2u, 2u, c), std::sqrt(15.0 / (32 * PI)) * s * s, 1e-15), "Y22(%g)", theta);
        }

#if defined(__cpp_lib_math_special_functions)
        // std::sph_legendre(l, m, θ) = Yl,m(θ, 0), including the Condon-Shortley phase
        for (unsigned l = 0; l < QuantumState::max_n; ++l)
            for (unsigned m = 0; m <= l; ++m)
                for (real_t theta : { 0.0, 0.1, 0.7, 1.3, PI / 2, 2.2, 3.0, PI })
                {
                    const real_t expected = std::sph_legendre(l, m, theta);
                    check(close(ylm_assoc_legendre(l, m, std::cos(theta)), expected, 1e-12),
                          "P~(l=%u, m=%u, theta=%g) = %.17g, std::sph_legendre = %.17g", l, m, theta,
                          ylm_assoc_legendre(l, m, std::cos(theta)), expected);
                }
#else
        std::printf("  (std::sph_legendre not available: compared with the closed forms only)\n");
#endif
    }

    ///////////////////////////////////////////////////////////////////
    ///
    ///     spherical_harmonic, spherical_real_basis: orthonormality
    ///

    void test_spherical_harmonics()
    {
        constexpr int L = 6;

        // θ: Clenshaw–Curtis, exact for polynomials in cosθ of degree <= N (the products here have degree <= 2L);
        // ϕ: periodic trapezoid rule, exact for e^(ikϕ) with |k| < Nphi (here |k| <= 2L).
        constexpr std::size_t N = 2 * L + 2;
        constexpr std::size_t Nphi = 4 * L;
        const auto w = Mesh::clenshaw_curtis_weights(N);

        auto inner = [&](auto&& f, auto&& g)
        {
            std::complex<real_t> sum{};
            for (std::size_t j = 0; j <= N; ++j)
                for (std::size_t k = 0; k < Nphi; ++k)
                {
                    const real_t theta = PI * static_cast<real_t>(j) / N;
                    const real_t phi = 2 * PI * static_cast<real_t>(k) / Nphi;
                    sum += w[j] * (2 * PI / Nphi) * f(theta, phi) * std::conj(g(theta, phi));
                }
            return sum;
        };

        for (int l1 = 0; l1 <= L; ++l1)
            for (int m1 = -l1; m1 <= l1; ++m1)
                for (int l2 = 0; l2 <= L; ++l2)
                    for (int m2 = -l2; m2 <= l2; ++m2)
                    {
                        const real_t delta = (l1 == l2 && m1 == m2) ? 1.0 : 0.0;
                        const auto u1 = static_cast<unsigned>(l1), u2 = static_cast<unsigned>(l2);

                        const auto c = inner(
                            [&](real_t t, real_t p) { return spherical_harmonic(u1, m1, t, p); },
                            [&](real_t t, real_t p) { return spherical_harmonic(u2, m2, t, p); });
                        check(std::abs(c - delta) < 1e-13, "<Y(%d,%d), Y(%d,%d)> = %.3g%+.3gi",
                              l1, m1, l2, m2, c.real(), c.imag());

                        const auto r = inner(
                            [&](real_t t, real_t p) { return std::complex<real_t>(spherical_real_basis(u1, m1, t, p)); },
                            [&](real_t t, real_t p) { return std::complex<real_t>(spherical_real_basis(u2, m2, t, p)); });
                        check(std::abs(r - delta) < 1e-13, "<Yreal(%d,%d), Yreal(%d,%d)> = %.3g",
                              l1, m1, l2, m2, r.real());
                    }
    }

    ///////////////////////////////////////////////////////////////////
    ///
    ///     QuantumState::Radial, QuantumState::RadialExtent
    ///

    void test_radial_extent()
    {
        // ∫_0^Rmax Rn,l(r)² r² dr (composite Simpson rule, fine enough to be exact to ~1e-12)
        // must be 1 minus the tail beyond Rmax, which RadialExtent() promises to be < 5e-6.
        constexpr std::size_t M = 20000;
        real_t worst_tail = 0.0;
        int worst_n = 0, worst_l = 0;

        for (int n = 1; n <= QuantumState::max_n; ++n)
            for (int l = 0; l < n; ++l)
            {
                const QuantumState qs{ n, l, 0 };
                const real_t Rmax = qs.RadialExtent();
                const real_t h = Rmax / M;

                real_t sum = 0.0;
                for (std::size_t i = 0; i <= M; ++i)
                {
                    const real_t r = h * static_cast<real_t>(i);
                    const real_t f = qs.Radial(r) * qs.Radial(r) * r * r;
                    sum += (i == 0 || i == M) ? f : (i % 2 ? 4 * f : 2 * f);
                }
                const real_t inside = sum * h / 3;
                const real_t tail = 1.0 - inside;

                check(tail > -1e-10 && tail < 5e-6, "n=%d, l=%d: probability beyond RadialExtent() = %.3g", n, l, tail);
                if (tail > worst_tail) { worst_tail = tail; worst_n = n; worst_l = l; }
            }

        std::printf("  largest tail beyond RadialExtent(): %.3g (n=%d, l=%d)\n", worst_tail, worst_n, worst_l);
    }

    ///////////////////////////////////////////////////////////////////
    ///
    ///     Mesh::for_state, Mesh::total_probability, Mesh::mean_r, mean_r2, mean_inv_r
    ///

    /// All states up to n = 6 in both bases, and a few large ones (the largest meshes).
    std::vector<QuantumState> test_states()
    {
        std::vector<QuantumState> states;
        for (Basis basis : { Basis::Complex, Basis::Real })
        {
            for (int n = 1; n <= 6; ++n)
                for (int l = 0; l < n; ++l)
                    for (int m = -l; m <= l; ++m)
                        states.emplace_back(n, l, m, basis);
            states.emplace_back(20, 0, 0, basis);
            states.emplace_back(20, 10, -5, basis);
            states.emplace_back(20, 19, 19, basis);
            states.emplace_back(20, 19, -19, basis);
        }
        return states;
    }

    const char* basis_name(const QuantumState& qs) { return qs.basis == Basis::Real ? "real" : "complex"; }

    void test_mesh_probability()
    {
        // The deficit is the radial quadrature error (largest for n = 1: ~3e-5) plus the tail beyond the mesh (< 5e-6);
        // the angular integrals are exact in both bases (for the real basis this needs > 2|m| azimuthal intervals)
        real_t worst = 0.0;
        const auto states = test_states();
        for (const auto& qs : states)
        {
            const Mesh mesh = Mesh::for_state(qs);
            const real_t deficit = std::abs(1.0 - mesh.total_probability);
            check(deficit < 5e-5, "n=%u, l=%u, m=%d (%s): total probability = %.12f",
                  qs.n, qs.l, qs.m, basis_name(qs), mesh.total_probability);
            worst = std::max(worst, deficit);
        }
        std::printf("  %zu states, largest |1 - P| = %.3g\n", states.size(), worst);
    }

    void test_expectation_values()
    {
        // <r>, <r²>, <1/r> (and the energy E = -<1/r>/2) on the default mesh agree with the exact values to ~1e-4
        // for n = 1 (radial step 0.2: 5 nodes per decay length of 1s), and better for larger n.
        // (With the plain trapezoid rule, <1/r> of 1s was 1.3% off; see Mesh::gregory_weights.)
        real_t worst_r = 0.0, worst_r2 = 0.0, worst_inv_r = 0.0;
        for (const auto& qs : test_states())
        {
            const Mesh mesh = Mesh::for_state(qs);
            const real_t err_r = std::abs(mesh.mean_r / qs.MeanR() - 1.0);
            const real_t err_r2 = std::abs(mesh.mean_r2 / qs.MeanR2() - 1.0);
            const real_t err_inv_r = std::abs(mesh.mean_inv_r / qs.MeanInvR() - 1.0);

            check(err_r < 1e-4, "n=%u, l=%u, m=%d (%s): <r> = %.12g, exact %.12g", qs.n, qs.l, qs.m, basis_name(qs), mesh.mean_r, qs.MeanR());
            check(err_r2 < 1e-4, "n=%u, l=%u, m=%d (%s): <r^2> = %.12g, exact %.12g", qs.n, qs.l, qs.m, basis_name(qs), mesh.mean_r2, qs.MeanR2());
            check(err_inv_r < 1e-4, "n=%u, l=%u, m=%d (%s): <1/r> = %.12g, exact %.12g", qs.n, qs.l, qs.m, basis_name(qs), mesh.mean_inv_r, qs.MeanInvR());

            worst_r = std::max(worst_r, err_r);
            worst_r2 = std::max(worst_r2, err_r2);
            worst_inv_r = std::max(worst_inv_r, err_inv_r);
        }
        std::printf("  largest relative errors: <r> %.3g, <r^2> %.3g, <1/r> %.3g\n", worst_r, worst_r2, worst_inv_r);

        // The exact values themselves: E = -<1/r>/2 (virial theorem), and two textbook cases
        const QuantumState s1{ 1, 0, 0 }, p2{ 2, 1, 0 };
        check(close(s1.MeanR(), 1.5, 1e-15) && close(s1.MeanR2(), 3.0, 1e-15) && close(s1.Energy(), -0.5, 1e-15), "1s exact values");
        check(close(p2.MeanR(), 5.0, 1e-15) && close(p2.MeanR2(), 30.0, 1e-15) && close(p2.Energy(), -0.125, 1e-15), "2p exact values");
        for (const auto& qs : { s1, p2, QuantumState{ 7, 3, 2 } })
            check(close(qs.Energy(), -qs.MeanInvR() / 2, 1e-15), "n=%u: E != -<1/r>/2", qs.n);
    }

    void test_gregory_weights()
    {
        // With p end corrections the rule is exact for polynomials of degree <= p (here on [0, 1], N intervals);
        // for few intervals the number of corrections drops to N/2, so the exactness to degree N/2.
        for (std::size_t N : { std::size_t{ 2 }, std::size_t{ 5 }, std::size_t{ 12 }, std::size_t{ 13 }, std::size_t{ 40 }, std::size_t{ 101 } })
        {
            const auto w = Mesh::gregory_weights(N);
            const real_t h = 1.0 / static_cast<real_t>(N);
            const auto p = static_cast<int>(std::min(Mesh::gregory_order, N / 2));
            for (int degree = 0; degree <= p; ++degree)
            {
                real_t sum = 0.0;
                for (std::size_t j = 0; j <= N; ++j)
                    sum += w[j] * h * std::pow(static_cast<real_t>(j) * h, degree);
                check(close(sum, 1.0 / (degree + 1), 1e-13), "N=%zu: integral of x^%d = %.17g, exact %.17g", N, degree, sum, 1.0 / (degree + 1));
            }
            check(*std::min_element(w.begin(), w.end()) > 0.0, "N=%zu: a non-positive weight", N);
        }
    }

    void test_radial_convergence()
    {
        // The Gregory rule (p = 6) converges much faster than the trapezoid rule (h²): halving the radial step divides the
        // error of the 1s integral ∫_0^R r·R(r)² dr (the <1/r> integrand, which rises linearly from r = 0) by ~150-250
        // (approaching 2^8 = 256 as h decreases; the trapezoid rule: 4).
        // Compared with the exact integral over the same range [0, R]: 1 - e^(-2R)·(1 + 2R).
        const QuantumState s1{ 1, 0, 0 };
        const real_t R = s1.RadialExtent();
        const real_t exact = 1.0 - std::exp(-2.0 * R) * (1.0 + 2.0 * R);

        real_t previous = 0.0;
        for (real_t h : { 0.2, 0.1, 0.05, 0.025 })
        {
            const Mesh mesh = Mesh::for_state(s1, Resolution{ h, 0, 0 });
            const real_t error = std::abs(mesh.mean_inv_r * mesh.total_probability - exact);
            std::printf("  h = %-6g <1/r> integral error %.3e%s\n", h, error,
                        previous > 0.0 ? ("  (" + std::to_string(previous / error) + "x smaller)").c_str() : "");
            if (previous > 1e-12)
                check(previous / error > 64.0, "h = %g: the error fell only %.1fx (expected > 64x)", h, previous / error);
            previous = error;
        }
        check(previous < 1e-10, "h = 0.025: error %.3g (expected < 1e-10)", previous);
    }

    void test_resolution()
    {
        const QuantumState qs{ 3, 2, -2 };

        // A custom resolution gives the requested mesh
        const Mesh mesh = Mesh::for_state(qs, Resolution{ 0.5, 31, 61 });
        check(mesh.psi.rows() == 31 && mesh.psi.columns() == 61, "custom angular points: %zu x %zu", mesh.psi.rows(), mesh.psi.columns());
        check(close(mesh.radial.delta, 0.5, 1e-12) && close(mesh.radial.back(), qs.RadialExtent(), 1e-15), "custom radial step");
        check(std::abs(1.0 - mesh.total_probability) < 1e-3, "custom mesh: total probability %.12f", mesh.total_probability);

        // Invalid resolutions are rejected
        auto rejected = [&](const Resolution& res)
        {
            try { Mesh::for_state(qs, res); }
            catch (const std::invalid_argument&) { return true; }
            return false;
        };
        check(rejected({ 0.0, 0, 0 }), "radial step 0 accepted");
        check(rejected({ -0.1, 0, 0 }), "negative radial step accepted");
        check(rejected({ 0.2, 20, 0 }), "even number of polar points accepted");
        check(rejected({ 0.2, 3, 0 }), "polar points < 2l + 1 accepted");
        check(rejected({ 0.2, 0, 40 }), "odd number of azimuthal intervals accepted");
        check(rejected({ 0.2, 0, 5 }), "azimuthal intervals <= 2|m| accepted");
        check(rejected({ 0.001, 1001, 2001 }), "a mesh larger than Resolution::max_psi_bytes accepted");
        check(!rejected({ 0.2, 5, 7 }), "the smallest valid angular mesh for l = 2, |m| = 2 rejected");
    }

    void test_orbital_names()
    {
        struct { int n, l, m; Basis basis; const char* name; } cases[] = {
            { 1, 0, 0, Basis::Complex, "1s" },        { 1, 0, 0, Basis::Real, "1s" },
            { 2, 1, 1, Basis::Real, "2p_x" },         { 2, 1, -1, Basis::Real, "2p_y" },     { 2, 1, 0, Basis::Real, "2p_z" },
            { 3, 2, -2, Basis::Real, "3d_xy" },       { 3, 2, 0, Basis::Real, "3d_z2" },     { 3, 2, 2, Basis::Real, "3d_x2-y2" },
            { 4, 3, -2, Basis::Real, "4f_xyz" },      { 4, 3, 3, Basis::Real, "4f_x(x2-3y2)" },
            { 5, 4, -2, Basis::Real, "5g_m-2" },      { 4, 3, 1, Basis::Complex, "4f(m=1)" },
            { 20, 19, -19, Basis::Complex, "20y(m=-19)" },
        };
        for (const auto& c : cases)
        {
            const auto name = QuantumState{ c.n, c.l, c.m, c.basis }.Name();
            check(name == c.name, "n=%d, l=%d, m=%d: name %s, expected %s", c.n, c.l, c.m, name.c_str(), c.name);
        }
    }
}

int main(int argc, char** argv)
{
    const std::map<std::string, std::function<void()>> tests
    {
        { "laguerre", test_laguerre },
        { "legendre", test_legendre },
        { "spherical_harmonics", test_spherical_harmonics },
        { "radial_extent", test_radial_extent },
        { "mesh_probability", test_mesh_probability },
        { "expectation_values", test_expectation_values },
        { "gregory_weights", test_gregory_weights },
        { "radial_convergence", test_radial_convergence },
        { "resolution", test_resolution },
        { "orbital_names", test_orbital_names },
    };

    std::vector<std::string> names(argv + 1, argv + argc);
    if (names.empty())
        for (const auto& [name, test] : tests) names.push_back(name);

    for (const auto& name : names)
    {
        const auto it = tests.find(name);
        if (it == tests.end())
        {
            std::printf("Unknown test: %s\n", name.c_str());
            return EXIT_FAILURE;
        }
        const int before = g_failures;
        std::printf("%s\n", name.c_str());
        it->second();
        std::printf("  %s\n", g_failures == before ? "passed" : "FAILED");
    }

    return g_failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
