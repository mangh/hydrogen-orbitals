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
    ///     Mesh::for_state, Mesh::compute_total_probability
    ///

    void test_mesh_probability()
    {
        // All states up to n = 6, and a few large ones (the largest meshes).
        std::vector<QuantumState> states;
        for (int n = 1; n <= 6; ++n)
            for (int l = 0; l < n; ++l)
                for (int m = -l; m <= l; ++m)
                    states.emplace_back(n, l, m);
        states.emplace_back(20, 0, 0);
        states.emplace_back(20, 10, -5);
        states.emplace_back(20, 19, 19);

        // The deficit is the radial trapezoid error (largest for n = 1: ~1e-4) plus the tail beyond the mesh (< 5e-6)
        real_t worst = 0.0;
        for (const auto& qs : states)
        {
            const Mesh mesh = Mesh::for_state(qs);
            const real_t deficit = std::abs(1.0 - mesh.total_probability);
            check(deficit < 2e-4, "n=%u, l=%u, m=%d: total probability = %.12f", qs.n, qs.l, qs.m, mesh.total_probability);
            worst = std::max(worst, deficit);
        }
        std::printf("  %zu states, largest |1 - P| = %.3g\n", states.size(), worst);
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
