#ifndef MANGH_HYDROGEN_SOLUTION_DATA_HPP
#define MANGH_HYDROGEN_SOLUTION_DATA_HPP

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <ostream>
#include <vector>

#include "math/constants.hpp"
#include "linspace.hpp"
#include "numbers.hpp"
#include "npy.hpp"
#include "array3d.hpp"
#include "quantum_state.hpp"

#define KahanSummation(total, addend, error)        \
                real_t addend_ = addend - error;    \
                real_t total_ = total + addend_;    \
                error = (total_ - total) - addend_; \
                total = total_;

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
        ///     Constants
        ///

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

        const QuantumState qs;

        /// @brief Total probability (or normalizing factor) for the entire area covered by the coordinates (r, θ, φ).
        real_t normalizing_factor;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        Mesh() = delete;

        /**
         * @brief Mesh constructor.
         *
         *
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
            normalizing_factor = compute_total_probability();
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        /// @brief Computes the wave function (Ψ) values on the mesh.
        void compute_wave_function()
        {
            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = 0; phi < azimuthal.size(); ++phi)
                    {
                        psi[r][theta][phi] = qs.WaveFunction(radial[r], polar[theta], azimuthal[phi]);
                    }
                }
            }
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
         * so the polar integral is exact for every l <= N/2 (l <= 9 for the 20-point polar grid).
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
                    s += b / (4.0 * k * k - 1.0) * std::cos(2.0 * k * theta);
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
        real_t compute_total_probability()
        {
            const auto w_theta = clenshaw_curtis_weights(polar.size() - 1);

            real_t total_probability{};

            // Variable to store the error (Kahan summation algorithm)
            real_t error = 0.0;

            for (std::size_t r = 1; r < radial.size(); ++r)
            {
                const real_t w_r = radial[r] * radial[r] * radial.delta * (r + 1 == radial.size() ? 0.5 : 1.0);

                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    const real_t w_rt = w_r * w_theta[theta] * azimuthal.delta;

                    for (std::size_t phi = 0; phi + 1 < azimuthal.size(); ++phi)
                    {
                        KahanSummation(total_probability, w_rt * std::norm(psi[r][theta][phi]), error)
                    }
                }
            }

            return total_probability;
        }

        /**
         * @brief Saves the mesh as NumPy (.npy) files in the directory @p dir:
         *
         * - `r.npy`, `theta.npy`, `phi.npy` - coordinates (r, θ, ϕ),
         * - `psi.npy` - probability amplitude Ψ, complex array of shape (r, θ, ϕ),
         * - `qn.npy`  - quantum numbers [n, l, m].
         *
         * @param dir - output directory (must exist).
         */
        void save(const std::filesystem::path& dir) const
        {
            npy::save(dir / "r.npy", radial.data(), { radial.size() });
            npy::save(dir / "theta.npy", polar.data(), { polar.size() });
            npy::save(dir / "phi.npy", azimuthal.data(), { azimuthal.size() });

            npy::save<complex_t>(dir / "psi.npy", { radial.size(), polar.size(), azimuthal.size() },
                [&](std::ostream& os)
                {
                    for (std::size_t r = 0; r < radial.size(); ++r)
                        for (std::size_t theta = 0; theta < polar.size(); ++theta)
                            npy::write_data(os, psi[r][theta].data(), azimuthal.size());
                });

            const std::int32_t qn[] = { static_cast<std::int32_t>(qs.n), static_cast<std::int32_t>(qs.l), qs.m };
            npy::save(dir / "qn.npy", qn, { 3 });
        }
    };
}}

#endif /* !MANGH_HYDROGEN_SOLUTION_DATA_HPP */
