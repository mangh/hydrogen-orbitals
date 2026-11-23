#ifndef MANGH_HYDROGENIUM_MESH_HPP
#define MANGH_HYDROGENIUM_MESH_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "math/constants.hpp"
#include "math/kahan_sum.hpp"
#include "linspace.hpp"
#include "numbers.hpp"
#include "npy.hpp"
#include "array3d.hpp"
#include "quantum_state.hpp"

namespace mangh { namespace hydrogenium
{
    using namespace mangh::math;

    /**
     * @brief Mesh resolution; a zero means the default (see Mesh::for_state).
     */
    struct Resolution
    {
        /// @brief Radial step (in a0).
        real_t radial_step = 0.2;

        /// @brief Number of polar nodes θ ∈ [0, π] (odd, >= 2l + 1); 0: max(21, 2l + 1).
        std::size_t polar_points = 0;

        /// @brief Number of azimuthal nodes ϕ ∈ [0, 2π], the duplicate node 2π included
        /// (the number of intervals, points - 1, must be even, >= 4 and > 2|m|); 0: max(40, 4|m|) + 1.
        std::size_t azimuthal_points = 0;

        /// @brief The largest mesh allowed: the size of psi in bytes (4 GB).
        static constexpr std::size_t max_psi_bytes = std::size_t{4000} * 1000 * 1000;
    };

    /**
     * @brief Data grid (in spherical coordinates)
     *
     */
    struct Mesh
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Radial coordinate.
        linspace<real_t> radial;

        /// @brief Polar coordinate θ (theta): 0 ≤ θ ≤ π (angle from the positive Z axis).
        linspace<real_t> polar;

        /// @brief Azimuthal coordinate ϕ (phi): 0 ≤ ϕ < 2π (angle from the X-axis in the XY-plane).
        linspace<real_t> azimuthal;

        /// @brief Probability amplitude (the codomain of the wave function Ψ in spherical coordinates).
        array3d<complex_t> psi;

        /// @brief Quantum state (n, l, m) the mesh is computed for.
        const QuantumState qs;

        /// @brief Total probability ∫|Ψ|² dV over the volume covered by the coordinates (r, θ, φ).
        real_t total_probability;

        /// @brief Expectation values `<r>`, `<r²>`, `<1/r>` integrated over the mesh (and divided by total_probability),
        /// to be compared with the exact QuantumState::MeanR(), MeanR2(), MeanInvR().
        real_t mean_r;
        real_t mean_r2;
        real_t mean_inv_r;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        Mesh() = delete;

        /**
         * @brief Mesh constructor: computes Ψ of the quantum state @p qs on the mesh and integrates |Ψ|².
         *
         * @param radial    - radial nodes r (evenly spaced, starting at 0),
         * @param polar     - polar nodes θ (evenly spaced on [0, π], both poles included: see clenshaw_curtis_weights),
         * @param azimuthal - azimuthal nodes ϕ (evenly spaced on [0, 2π], the last node duplicating the first),
         * @param qs        - quantum state (n, l, m).
         */
        Mesh(const linspace<real_t>& radial,
             const linspace<real_t>& polar,
             const linspace<real_t>& azimuthal,
             const QuantumState& qs
            )
            : radial(radial)
            , polar(polar)
            , azimuthal(azimuthal)
            , psi(radial.size(), polar.size(), azimuthal.size())
            , qs(qs)
        {
            compute_wave_function();

            // The following calculations aren't necessary for the visualization;
            // they check the wave function (normalization, expectation values) and the quadrature:
            const auto density = radial_density();
            total_probability = integrate_radial(density, [](real_t) { return 1.0; });
            mean_r = integrate_radial(density, [](real_t r) { return r; }) / total_probability;
            mean_r2 = integrate_radial(density, [](real_t r) { return r * r; }) / total_probability;
            mean_inv_r = integrate_radial(density, [](real_t r) { return 1.0 / r; }) / total_probability;
        }

        /**
         * @brief The mesh used by the application for the quantum state @p qs.
         *
         * - r ∈ [0, RadialExtent()] (the range scales with n, see QuantumState::RadialExtent), with the step @p res.radial_step
         *   (rounded so that it divides the range evenly),
         *
         * - θ ∈ [0, π]: an odd number of nodes, so that a node lies exactly at θ = π/2 (XY-plane plots),
         *   and N = nodes - 1 >= 2l intervals, so that the θ integral (Clenshaw–Curtis, see clenshaw_curtis_weights) is exact for |Yl,m|²;
         *   by default max(21, 2l + 1) nodes,
         *
         * - ϕ ∈ [0, 2π] (the last node duplicates the first): an even number of intervals, so that for every node ϕ
         *   there is also a node at ϕ + π (the opposite half-plane in the XZ-plane plots), and > 2|m| intervals, which keeps
         *   the ϕ integral exact for real orbitals, where |Ψ|² contains cos²(mϕ); by default max(40, 4|m|) intervals,
         *   i.e. at least 4 samples per period of cos(mϕ) in the XY-plane plots.
         *
         * @exception std::invalid_argument - when @p res does not meet these conditions.
         */
        static Mesh for_state(const QuantumState& qs, const Resolution& res = {})
        {
            static constexpr real_t PI = mangh::math::PI_V<real_t>;

            auto str = [](real_t x) { std::ostringstream s; s << x; return s.str(); };

            const real_t radialMax = qs.RadialExtent();
            if (!(res.radial_step > 0.0 && res.radial_step <= radialMax / 2))
                throw std::invalid_argument("radial step = " + str(res.radial_step)
                    + " must satisfy 0 < step <= " + str(radialMax / 2) + " (half the radial range)");
            const auto radialPoints = static_cast<std::size_t>(std::lround(radialMax / res.radial_step)) + 1;

            const std::size_t polarMin = std::max<std::size_t>(3, 2 * qs.l + 1);
            const std::size_t polarPoints = res.polar_points ? res.polar_points : std::max<std::size_t>(21, polarMin);
            if (polarPoints % 2 == 0 || polarPoints < polarMin)
                throw std::invalid_argument("polar points = " + std::to_string(polarPoints)
                    + " must be odd and >= " + std::to_string(polarMin) + " (max(3, 2l + 1))");

            const auto twoM = static_cast<std::size_t>(2 * std::abs(qs.m));
            const std::size_t azimuthalPoints = res.azimuthal_points ? res.azimuthal_points : std::max<std::size_t>(40, 2 * twoM) + 1;
            const std::size_t azimuthalIntervals = azimuthalPoints > 0 ? azimuthalPoints - 1 : 0;
            if (azimuthalIntervals % 2 != 0 || azimuthalIntervals < 4 || azimuthalIntervals <= twoM)
                throw std::invalid_argument("azimuthal points = " + std::to_string(azimuthalPoints)
                    + " must be odd (an even number of intervals), >= 5 and > 2|m| + 1 = " + std::to_string(twoM + 1));

            const std::size_t bytes = radialPoints * polarPoints * azimuthalPoints * sizeof(complex_t);
            if (bytes > Resolution::max_psi_bytes)
                throw std::invalid_argument("the mesh " + std::to_string(radialPoints) + " x " + std::to_string(polarPoints)
                    + " x " + std::to_string(azimuthalPoints) + " would need " + std::to_string(bytes / 1000000)
                    + " MB for psi (limit: " + std::to_string(Resolution::max_psi_bytes / 1000000)
                    + " MB): use a larger radial step or fewer angular points");

            return Mesh
            {
                linspace<real_t>{ 0.0, radialMax, radialPoints },     //  radial coordinates (r)
                linspace<real_t>{ 0.0, PI, polarPoints },             //  polar coordinates (θ)
                linspace<real_t>{ 0.0, 2.0 * PI, azimuthalPoints },   //  azimuthal coordinates (ϕ)
                qs
            };
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        /**
         * @brief Computes the wave function (Ψ) values on the mesh.
         *
         * @remarks Ψ(r,θ,ϕ) = Rn,l(r) · Yl,m(θ,ϕ) is separable, so Rn,l is evaluated once per radial node
         * and Yl,m once per angular node (θ,ϕ); each mesh value is then a single product
         * (the same product as QuantumState::WaveFunction, hence the same result).
         */
        void compute_wave_function()
        {
            std::vector<real_t> R(radial.size());
            for (std::size_t r = 0; r < radial.size(); ++r)
                R[r] = qs.Radial(radial[r]);

            std::vector<complex_t> Y(polar.size() * azimuthal.size());
            for (std::size_t theta = 0; theta < polar.size(); ++theta)
                for (std::size_t phi = 0; phi < azimuthal.size(); ++phi)
                    Y[theta * azimuthal.size() + phi] = qs.Angular(polar[theta], azimuthal[phi]);

            for (std::size_t r = 0; r < radial.size(); ++r)
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                    for (std::size_t phi = 0; phi < azimuthal.size(); ++phi)
                        psi(r, theta, phi) = R[r] * Y[theta * azimuthal.size() + phi];
        }

        /**
         * @brief Clenshaw–Curtis quadrature weights for the polar integral.
         *
         * - ∫_0^π f(θ) sinθ dθ  =  ∫_-1^1 f(x) dx  ≈  Σ w_j f(θ_j),   x = cosθ,   θ_j = jπ/N,  j = 0..N.
         *
         * @param N - number of polar intervals (= polar.size() - 1); the nodes must be evenly spaced on [0, π], poles included.
         *
         * @returns N + 1 weights.
         *
         * @remarks The rule interpolates f(x) at the nodes x_j = cos(jπ/N) by a polynomial of degree N
         * (equivalently, a cosine series in θ) and integrates that polynomial exactly. Hence the result is exact
         * whenever f is a polynomial in cosθ of degree <= N. |Yl,m(θ,ϕ)|² is such a polynomial of degree 2l,
         * so the polar integral is exact for every l <= N/2 (l <= 10 for the default 21-node polar grid).
         *
         * Weights (Trefethen, "Is Gauss quadrature better than Clenshaw–Curtis?", SIAM Review 50 (2008)):
         *
         * - w_j = (c_j / N) · (1 − Σ_{k=1}^{⌊N/2⌋} b_k / (4k² − 1) · cos(2kθ_j)),
         *
         * - c_0 = c_N = 1, otherwise c_j = 2;  b_k = 1 for k = N/2, otherwise b_k = 2.
         *
         */
        static std::vector<real_t> clenshaw_curtis_weights(std::size_t N)
        {
            static constexpr real_t PI = mangh::math::PI_V<real_t>;

            std::vector<real_t> w(N + 1);
            for (std::size_t j = 0; j <= N; ++j)
            {
                const real_t theta = PI * static_cast<real_t>(j) / static_cast<real_t>(N);

                real_t s = 0.0;
                for (std::size_t k = 1; k <= N / 2; ++k)
                {
                    const real_t b = (2 * k == N) ? 1.0 : 2.0;
                    const auto kr = static_cast<real_t>(k);
                    s += b / (4.0 * kr * kr - 1.0) * std::cos(2.0 * kr * theta);
                }

                const real_t c = (j == 0 || j == N) ? 1.0 : 2.0;
                w[j] = c / static_cast<real_t>(N) * (1.0 - s);
            }
            return w;
        }

        /**
         * @brief Integrates |Ψ|² over the angles at every radial node: `ρ(r) = ∫ |Ψ(r,θ,ϕ)|² sinθ dθ dϕ`.
         *
         * - θ: Clenshaw–Curtis rule (exact for |Yl,m|² with l <= (polar.size() - 1) / 2),
         * - ϕ: periodic trapezoid rule (the duplicate node ϕ = 2π is skipped).
         *
         * @returns ρ at the radial nodes (for an exact quadrature, ρ(r) = Rn,l(r)²).
         */
        std::vector<real_t> radial_density() const
        {
            const auto w_theta = clenshaw_curtis_weights(polar.size() - 1);

            std::vector<real_t> density(radial.size());
            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                // Kahan summation keeps the rounding error of the many small terms in check
                KahanSum<real_t> sum;
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                    for (std::size_t phi = 0; phi + 1 < azimuthal.size(); ++phi)
                        sum += w_theta[theta] * azimuthal.delta * std::norm(psi(r, theta, phi));
                density[r] = sum.value();
            }
            return density;
        }

        /**
         * @brief Radial integral `∫ f(r) ρ(r) r² dr` over the mesh (trapezoid rule), i.e. `∫ f(r) |Ψ|² dV`.
         *
         * @param density - ρ(r) at the radial nodes (see radial_density),
         * @param f       - function of r; it is not evaluated at r = 0, where r² makes the integrand vanish
         *                  for f = 1, r, r², 1/r (but not for 1/r²).
         *
         * @returns The integral (for f = 1: the total probability, which should be very close to 1 for a normalized Ψ).
         */
        template<typename F>
        real_t integrate_radial(const std::vector<real_t>& density, F&& f) const
        {
            KahanSum<real_t> total;
            for (std::size_t r = 1; r < radial.size(); ++r)
            {
                const real_t w_r = radial[r] * radial[r] * radial.delta * (r + 1 == radial.size() ? 0.5 : 1.0);
                total += w_r * f(radial[r]) * density[r];
            }
            return total.value();
        }

        /**
         * @brief Saves the mesh as NumPy (.npy) files in the directory @p dir:
         *
         * - `r.npy`, `theta.npy`, `phi.npy` - coordinates (r, θ, ϕ),
         * - `psi.npy` - probability amplitude Ψ, complex array of shape (r, θ, ϕ),
         * - `qn.npy`  - quantum numbers [n, l, m],
         * - `info.json` - the state (n, l, m, basis, orbital name), the mesh shape, and the checks
         *   (total probability; expectation values and energy, on the mesh and exact).
         *
         * @param dir - output directory (must exist).
         *
         * @exception std::runtime_error - when a file cannot be written.
         */
        void save(const std::filesystem::path& dir) const
        {
            npy::save(dir / "r.npy", radial.data(), { radial.size() });
            npy::save(dir / "theta.npy", polar.data(), { polar.size() });
            npy::save(dir / "phi.npy", azimuthal.data(), { azimuthal.size() });

            npy::save(dir / "psi.npy", psi.data(), { psi.pages(), psi.rows(), psi.columns() });

            const std::int32_t qn[] = { static_cast<std::int32_t>(qs.n), static_cast<std::int32_t>(qs.l), qs.m };
            npy::save(dir / "qn.npy", qn, { 3 });

            const auto file = dir / "info.json";
            std::ofstream os(file);
            if (!os) throw std::runtime_error("Cannot open " + file.string() + " for writing");
            os << info_json();
            os.close();
            if (!os) throw std::runtime_error("Failed to write " + file.string());
        }

        /// @brief The contents of info.json (see save); the energy on the mesh comes from the virial theorem, E = -<1/r>/2.
        std::string info_json() const
        {
            auto number = [](real_t x)
            {
                std::ostringstream s;
                s << std::setprecision(std::numeric_limits<real_t>::max_digits10) << x;
                return s.str();
            };
            auto pair = [&](real_t grid, real_t exact)
            {
                return "{ \"grid\": " + number(grid) + ", \"exact\": " + number(exact) + " }";
            };

            std::ostringstream js;
            js << "{\n"
               << "  \"n\": " << qs.n << ",\n"
               << "  \"l\": " << qs.l << ",\n"
               << "  \"m\": " << qs.m << ",\n"
               << "  \"basis\": \"" << (qs.basis == Basis::Real ? "real" : "complex") << "\",\n"
               << "  \"orbital\": \"" << qs.Name() << "\",\n"
               << "  \"shape\": [" << psi.pages() << ", " << psi.rows() << ", " << psi.columns() << "],\n"
               << "  \"radial_step\": " << number(radial.delta) << ",\n"
               << "  \"radial_max\": " << number(radial.back()) << ",\n"
               << "  \"total_probability\": " << number(total_probability) << ",\n"
               << "  \"mean_r\": " << pair(mean_r, qs.MeanR()) << ",\n"
               << "  \"mean_r2\": " << pair(mean_r2, qs.MeanR2()) << ",\n"
               << "  \"mean_inv_r\": " << pair(mean_inv_r, qs.MeanInvR()) << ",\n"
               << "  \"energy_hartree\": " << pair(-mean_inv_r / 2.0, qs.Energy()) << "\n"
               << "}\n";
            return js.str();
        }
    };
}}

#endif /* !MANGH_HYDROGENIUM_MESH_HPP */
