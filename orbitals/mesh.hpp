#ifndef MANGH_HYDROGEN_SOLUTION_DATA_HPP
#define MANGH_HYDROGEN_SOLUTION_DATA_HPP

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <ostream>
#include <vector>

#include "math/constants.hpp"
#include "linspace.hpp"
#include "numbers.hpp"
#include "array3d.hpp"
#include "quantum_state.hpp"

#define KahanSummation(total, addend, error)        \
                real_t addend_ = addend - error;  \
                real_t total_ = total + addend_;  \
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

        /// @brief Maximum probability density value.
        real_t max_pd;

        /// @brief Maximum value of the real part of the probability amplitude (Ψ).
        real_t max_real;

        /// @brief Total probability (or normalizing factor) for the entire area covered by the coordinates (r, θ, φ).
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
         * @brief Mesh constructor.
         *
         *
         */
        Mesh(const linspace<real_t>& radial,
             const linspace<real_t>& polar,
             const linspace<real_t>& azimuthal
            )
            : radial(radial)
            , polar(polar)
            , azimuthal(azimuthal)

            , psi(radial.size(), polar.size(), azimuthal.size())

            , max_pd(std::numeric_limits<real_t>::lowest())
            , max_real(std::numeric_limits<real_t>::lowest())
            , total_probability(std::numeric_limits<real_t>::lowest())
        {
            // // What's the point of the below code?
            // // The azimuthal coordinates are already evenly spaced in the interval [0, 2π) due to the linspace constructor.
            // // Adjusting them by half of the delta would shift the entire set of azimuthal coordinates by a constant value,
            // // which is not necessary for the purpose of defining a mesh in spherical coordinates.
            // // The azimuthal coordinate should start at 0 and end just before 2π to cover the full circle without overlap.
            // for(size_t phi = 0; phi < this->azimuthal.size(); ++phi)
            // {
            //     this->azimuthal[phi] += azimuthal.delta / 2.0;
            // }
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        /// @brief Cartesian coordinate (x) from spherical coordinates (r, θ, φ).
        real_t x(std::size_t r, std::size_t theta, std::size_t phi) const
        {
            return radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
        }
        
        /// @brief Cartesian coordinate (y) from spherical coordinates (r, θ, φ).
        real_t y(std::size_t r, std::size_t theta, std::size_t phi) const
        {
            return radial[r] * sin(polar[theta]) * sin(azimuthal[phi]);
        }

        /// @brief Cartesian coordinate (z) from spherical coordinates (r, θ, φ).
        real_t z(std::size_t r, std::size_t theta, [[maybe_unused]] std::size_t phi) const
        {
            return radial[r] * cos(polar[theta]);
        }

        void compute_wave_function(const QuantumState& qs)
        {
            max_pd = std::numeric_limits<real_t>::lowest();
            max_real = std::numeric_limits<real_t>::lowest();

            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = 0; phi < azimuthal.size(); ++phi)
                    {
                        // probability amplitude (Ψ) 
                        complex_t pa = qs.WaveFunction(radial[r], polar[theta], azimuthal[phi]);
                        if (pa.real() > max_real)
                            max_real = pa.real();

                        // probability density (|Ψ|^2)
                        real_t pd = std::norm(pa);
                        if (pd > max_pd)
                            max_pd = pd;

                        psi[r][theta][phi] = pa;
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
         * @brief Integrates |Ψ|² over the mesh (the result is stored in total_probability; ≈ 1 for a normalized Ψ).
         *
         * - r: trapezoid rule (r = 0 contributes nothing because of r²),
         * - θ: Clenshaw–Curtis rule (exact for |Yl,m|² with l <= (polar.size() - 1) / 2),
         * - ϕ: periodic trapezoid rule (the duplicate node ϕ = 2π is skipped).
         *
         */
        void normalize_psi()
        {
            const auto w_theta = clenshaw_curtis_weights(polar.size() - 1);

            total_probability = 0.0;

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

            ///* Normalize */
            //auto integral_sqrt = sqrt(total_probability);
            //for (std::size_t r = 1; r < radial.size(); ++r)
            //{
            //    for (std::size_t theta = 0; theta < polar.size(); ++theta)
            //    {
            //        for (std::size_t phi = 0; phi < azimuthal.size(); ++phi)
            //        {
            //            psi[r][theta][phi] /= integral_sqrt;
            //        }
            //    }
            //}
        }

        void export_grid_data(std::ostream &os) const
        {
            os
              << radial.size() << " "
              << polar.size() << " "
              << azimuthal.size() << " "
              << max_pd << " "
              << max_real
              << std::endl;
        }

        void export_probability_density(std::ostream &os) const
        {
            /* Export psi central nodes */
            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = 0; phi < 1; ++phi)
                    {
                        os << x(r, theta, phi) << " " << z(r, theta, phi) << " " << norm(psi[r][theta][phi]) << std::endl;
                    }
                }
            }

            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = azimuthal.size() / 2; phi < azimuthal.size() / 2 + 1; ++phi)
                    {
                        os << x(r, theta, phi) << " " << z(r, theta, phi) << " " << norm(psi[r][theta][phi]) << std::endl;
                    }
                }
            }
        }

        void export_psi_real(std::ostream &os) const
        {
            std::size_t phi_slice = 0; // Should be a value between 0 and azimuthal.size() / 2 - 2

            /* Export psi central nodes */
            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = phi_slice; phi < phi_slice + 1; ++phi)
                    {
                        os << x(r, theta, phi) << " " << z(r, theta, phi) << " " << psi[r][theta][phi].real() << std::endl;
                    }
                }
            }

            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = azimuthal.size() / 2 + phi_slice; phi < azimuthal.size() / 2 + phi_slice + 1; ++phi)
                    {
                        os << x(r, theta, phi) << " " << z(r, theta, phi) << " " << psi[r][theta][phi].real() << std::endl;
                    }
                }
            }
        }

        void export_psi_im(std::ostream &os) const
        {
            std::size_t phi_slice = 0; // Should be a value between 0 and polar.size() / 2 - 2

            /* Export psi central nodes */
            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = phi_slice; phi < phi_slice + 1; ++phi)
                    {
                        os << x(r, theta, phi) << " " << z(r, theta, phi) << " " << psi[r][theta][phi].imag() << std::endl;
                    }
                }
            }

            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = azimuthal.size() / 2 + phi_slice; phi < azimuthal.size() / 2 + phi_slice + 1; ++phi)
                    {
                        os << x(r, theta, phi) << " " << z(r, theta, phi) << " " << psi[r][theta][phi].imag() << std::endl;
                    }
                }
            }
        }

        void export_data_phi_2(std::ostream &os) const
        {
            /* Export psi central nodes */
            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = polar.size() / 2; theta < polar.size() / 2 + 1; ++theta)
                {
                    for (std::size_t phi = 0; phi < azimuthal.size(); ++phi)
                    {
                        os << x(r, theta, phi) << " " << y(r, theta, phi) << " " << psi[r][theta][phi].real() << std::endl;
                    }
                }
            }
        }

        void export_data_phi_2_im(std::ostream &os) const
        {
            /* Export psi central nodes */
            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = polar.size() / 2; theta < polar.size() / 2 + 1; ++theta)
                {
                    for (std::size_t phi = 0; phi < azimuthal.size(); ++phi)
                    {
                        os << x(r, theta, phi) << " " << y(r, theta, phi) << " " << psi[r][theta][phi].imag() << std::endl;
                    }
                }
            }
        }

    //    void verify_solution(const Parameters& params, const QuantumState& qs)
    //    {
    //        static constexpr complex_t CZERO{ 0.0, 0.0 };

    //        array3d<complex_t> ratio(radial.size(), azimuthal.size(), polar.size());

    //        real_t alpha = params.h / (2.0 * params.mp);
    //        real_t E = params.ke / (qs.n * qs.n);

    //        for (std::size_t r = 1; r < radial.size() - 1; ++r)
    //        {
    //            for (std::size_t theta = 1; theta < polar.size() - 1; ++theta)
    //            {
    //                for (std::size_t phi = 1; phi < azimuthal.size() - 1; ++phi)
    //                {
    //                    if (psi[r][theta][phi] == CZERO)
    //                        continue;

    //                    real_t dV1 = radial[r] * radial[r] * radial.delta;
    //                    real_t dV2 = radial[r] * radial[r] * sin(polar[theta]) * sin(polar[theta]) * azimuthal.delta;
    //                    real_t dV3 = radial[r] * radial[r] * sin(polar[theta]) * polar.delta;

    //                    // real_t aS = r[i - 1] * r[i - 1] * alpha / (dV1 * r_p.delta);
    //                    // real_t aN = r[i] * r[i] * alpha / (dV1 * r_p.delta);
    //                    real_t r1 = radial[r] - radial.delta / 2.0;
    //                    real_t r2 = radial[r] + radial.delta / 2.0;

    //                    real_t aS = r1 * r1 * alpha / (dV1 * radial.delta);
    //                    real_t aN = r2 * r2 * alpha / (dV1 * radial.delta);
    //                    real_t aW = alpha / (dV2 * azimuthal.delta);
    //                    real_t aE = alpha / (dV2 * azimuthal.delta);
    //                    real_t aB = sin(polar[theta] - polar.delta / 2.0) * alpha / (dV3 * polar.delta);
    //                    real_t aT = sin(polar[theta] + polar.delta / 2.0) * alpha / (dV3 * polar.delta);

    //                    complex_t psi_s{ -aS * psi[r - 1][phi][theta] };
    //                    complex_t psi_n{ -aN * psi[r + 1][phi][theta] };
    //                    complex_t psi_w{ -aW * psi[r][phi - 1][theta] };
    //                    complex_t psi_e{ -aE * psi[r][phi + 1][theta] };
    //                    complex_t psi_b{ -aB * psi[r][phi][theta - 1] };
    //                    complex_t psi_t{ -aT * psi[r][phi][theta + 1] };

    //                    complex_t psi_p_k_rp{ -params.kp / radial[r] * psi[r][phi][theta] };

    //                    complex_t Fp{ psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + psi_p_k_rp };
    //                    real_t ap = aS + aN + aW + aE + aB + aT - E;
    //                    complex_t ap_psi_p{ -ap * psi[r][theta][phi] };

    //                    ratio[r][theta][phi] = Fp / ap_psi_p;
    //                }
    //            }
    //        }

    //        /* Compute error */
    //        real_t error = 0.0;
    //        for (std::size_t ir = 1; ir < radial.size() - 1; ++ir)
    //        {
    //            for (std::size_t it = 1; it < polar.size() - 1; ++it)
    //            {
    //                for (std::size_t ip = 1; ip < azimuthal.size() - 1; ++ip)
    //                {
    //                    if (fabs(ratio[ir][it][ip].real()) > 0.2)
    //                    {
    //                        error += fabs(1.0 - ratio[ir][it][ip].real());
    //                    }
    //                }
    //            }
    //        }
    //        error = 100.0 * error / ((radial.size() - 2) * (azimuthal.size() - 2) * (polar.size() - 2));
    //    }
    };
}}

#endif /* !MANGH_HYDROGEN_SOLUTION_DATA_HPP */
