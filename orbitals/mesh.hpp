#ifndef MANGH_HYDROGENIUM_MESH_HPP
#define MANGH_HYDROGENIUM_MESH_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
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

            // The following calculation isn't necessary
            // when we use a normalized function Ψ (for which the total probability ≈ 1).
            // We do it just to check if that's really the case:
            total_probability = compute_total_probability();
        }

        /**
         * @brief The mesh used by the application for the quantum state @p qs.
         *
         * - r ∈ [0, RadialExtent()], with a fixed step of 0.2 (the range scales with n, see QuantumState::RadialExtent),
         *
         * - θ ∈ [0, π]: an odd number of nodes, so that a node lies exactly at θ = π/2 (XY-plane plots),
         *   and N = nodes - 1 >= 2l intervals, so that the θ integral (Clenshaw–Curtis, see compute_total_probability) is exact for |Yl,m|²,
         *
         * - ϕ ∈ [0, 2π] (the last node duplicates the first): an even number of intervals, so that for every node ϕ
         *   there is also a node at ϕ + π (the opposite half-plane in the XZ-plane plots), and >= 4|m| intervals:
         *   at least 4 samples per period of cos(mϕ) in the XY-plane plots (and > 2|m|, which also keeps the ϕ integral
         *   exact for real orbitals, where |Ψ|² contains cos²(mϕ)).
         */
        static Mesh for_state(const QuantumState& qs)
        {
            static constexpr real_t PI = mangh::math::PI_V<real_t>;

            constexpr real_t radialStep = 0.2;
            const real_t radialMax = qs.RadialExtent();
            const auto radialPoints = static_cast<std::size_t>(std::lround(radialMax / radialStep)) + 1;

            const std::size_t polarPoints = std::max<std::size_t>(21, 2 * qs.l + 1);

            const std::size_t azimuthalPoints = std::max<std::size_t>(40, static_cast<std::size_t>(4 * std::abs(qs.m))) + 1;

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
         * @brief Integrates |Ψ|² over the mesh.
         *
         * - r: trapezoid rule (r = 0 contributes nothing because of r²),
         * - θ: Clenshaw–Curtis rule (exact for |Yl,m|² with l <= (polar.size() - 1) / 2),
         * - ϕ: periodic trapezoid rule (the duplicate node ϕ = 2π is skipped).
         *
         * @remarks The result should be very close to 1 for a normalized Ψ.
         *
         * @returns The total probability for the whole mesh.
         */
        real_t compute_total_probability() const
        {
            const auto w_theta = clenshaw_curtis_weights(polar.size() - 1);

            // Kahan summation keeps the rounding error of the many small terms in check
            KahanSum<real_t> total;

            for (std::size_t r = 1; r < radial.size(); ++r)
            {
                const real_t w_r = radial[r] * radial[r] * radial.delta * (r + 1 == radial.size() ? 0.5 : 1.0);

                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    const real_t w_rt = w_r * w_theta[theta] * azimuthal.delta;

                    for (std::size_t phi = 0; phi + 1 < azimuthal.size(); ++phi)
                    {
                        total += w_rt * std::norm(psi(r, theta, phi));
                    }
                }
            }

            return total.value();
        }

        /**
         * @brief Saves the mesh as NumPy (.npy) files in the directory @p dir:
         *
         * - `r.npy`, `theta.npy`, `phi.npy` - coordinates (r, θ, ϕ),
         * - `psi.npy` - probability amplitude Ψ, complex array of shape (r, θ, ϕ),
         * - `qn.npy`  - quantum numbers [n, l, m].
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
        }
    };
}}

#endif /* !MANGH_HYDROGENIUM_MESH_HPP */
