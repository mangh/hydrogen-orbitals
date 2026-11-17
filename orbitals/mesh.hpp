#ifndef MANGH_HYDROGEN_SOLUTION_DATA_HPP
#define MANGH_HYDROGEN_SOLUTION_DATA_HPP

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <ostream>

#include "math/constants.hpp"
#include "linspace.hpp"
#include "numbers.hpp"
#include "array3d.hpp"
#include "quantum_state.hpp"

namespace mangh { namespace hydrogenium
{
    using namespace mangh::math;

    /**
     * @brief Data mesh (spherical coordinate grid)
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
        linspace<number_t> radial;

        /// @brief Polar coordinate θ (theta): 0 ≤ θ ≤ π (angle from the positive Z axis).
        linspace<number_t> polar;

        /// @brief Azimuthal coordinate ϕ (phi): 0 ≤ ϕ < 2π (angle from the X-axis in the XY-plane).
        linspace<number_t> azimuthal;

        /// @brief Wave function (Ψ).
        Array3d<complex_t> psi;

        /// @brief Maximum probability density value.
        number_t max_pd;

        /// @brief Maximum value of the real part of the wave function (Ψ).
        number_t max_real;

        /// @brief Total probability (or normalizing factor) for the entire area covered by the coordinates (r, θ, φ).
        number_t total_probability;

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
        Mesh(const linspace<number_t>& radial,
             const linspace<number_t>& polar,
             const linspace<number_t>& azimuthal
            )
            : radial(radial)
            , polar(polar)
            , azimuthal(azimuthal)

            , psi(radial.size(), polar.size(), azimuthal.size())

            , max_pd(std::numeric_limits<number_t>::lowest())
            , max_real(std::numeric_limits<number_t>::lowest())
            , total_probability(std::numeric_limits<number_t>::lowest())
        {
            for(auto& phi : this->azimuthal) phi += azimuthal.delta / 2.0;
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        void compute_wave_function(const QuantumState& qs)
        {
            max_pd = std::numeric_limits<number_t>::lowest();
            max_real = std::numeric_limits<number_t>::lowest();

            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = 0; phi < azimuthal.size(); ++phi)
                    {
                        complex_t c = qs.WaveFunction(radial[r], polar[theta], azimuthal[phi]);
                        if (c.real() > max_real)
                            max_real = c.real();

                        number_t pd = std::norm(c);
                        if (pd > max_pd)
                            max_pd = pd;

                        psi[r][theta][phi] = c;
                    }
                }
            }
        }

        void normalize_psi()
        {
            static constexpr number_t M_PI = mangh::math::PI_V<number_t>;

            total_probability = 0.0;

            /* Top pole */
            for (std::size_t r = 1; r < radial.size(); ++r)
            {
                number_t rlo = radial[r] - radial.delta / 2.0;
                number_t rhi = radial[r] + radial.delta / 2.0;
                number_t probability_per_vol =
                    2.0 / 3.0 * M_PI * (1 - cos(polar.delta / 2.0)) * (rhi * rhi * rhi - rlo * rlo * rlo) * norm(psi[r][0][0]);

                total_probability += probability_per_vol;
            }

            /* Central nodes */
            for (std::size_t r = 1; r < radial.size(); ++r)
            {
                for (std::size_t theta = 1; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = 0; phi < azimuthal.size() - 1; ++phi)
                    {
                        number_t probability_per_vol = 
                            radial[r] * radial[r] * sin(polar[theta]) * polar.delta * azimuthal.delta * radial.delta * norm(psi[r][theta][phi]);
                        total_probability += probability_per_vol;
                    }
                }
            }

            /* Bottom pole */
            for (std::size_t r = 1; r < radial.size(); ++r)
            {
                number_t rlo = radial[r] - radial.delta / 2.0;
                number_t rhi = radial[r] + radial.delta / 2.0;
                number_t probability_per_vol = 
                    2.0 / 3.0 * M_PI * (1 - cos(polar.delta / 2.0)) * (rhi * rhi * rhi - rlo * rlo * rlo) * norm(psi[r][polar.size() - 1][0]);
                total_probability += probability_per_vol;
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
                        number_t x = radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
                        number_t z = radial[r] * cos(polar[theta]);
                        os << x << " " << z << " " << norm(psi[r][theta][phi]) << std::endl;
                    }
                }
            }

            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = azimuthal.size() / 2; phi < azimuthal.size() / 2 + 1; ++phi)
                    {
                        number_t x = radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
                        number_t z = radial[r] * cos(polar[theta]);
                        os << x << " " << z << " " << norm(psi[r][theta][phi]) << std::endl;
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
                        number_t x = radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
                        number_t z = radial[r] * cos(polar[theta]);
                        os << x << " " << z << " " << psi[r][theta][phi].real() << std::endl;
                    }
                }
            }

            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = azimuthal.size() / 2 + phi_slice; phi < azimuthal.size() / 2 + phi_slice + 1; ++phi)
                    {
                        double x = radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
                        double z = radial[r] * cos(polar[theta]);
                        os << x << " " << z << " " << psi[r][theta][phi].real() << std::endl;
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
                        number_t x = radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
                        number_t z = radial[r] * cos(polar[theta]);
                        os << x << " " << z << " " << psi[r][theta][phi].imag() << std::endl;
                    }
                }
            }

            for (std::size_t r = 0; r < radial.size(); ++r)
            {
                for (std::size_t theta = 0; theta < polar.size(); ++theta)
                {
                    for (std::size_t phi = azimuthal.size() / 2 + phi_slice; phi < azimuthal.size() / 2 + phi_slice + 1; ++phi)
                    {
                        number_t x = radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
                        number_t z = radial[r] * cos(polar[theta]);
                        os << x << " " << z << " " << psi[r][theta][phi].imag() << std::endl;
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
                        number_t x = radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
                        number_t y = radial[r] * sin(polar[theta]) * sin(azimuthal[phi]);
                        os << x << " " << y << " " << psi[r][theta][phi].real() << std::endl;
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
                        number_t x = radial[r] * sin(polar[theta]) * cos(azimuthal[phi]);
                        number_t y = radial[r] * sin(polar[theta]) * sin(azimuthal[phi]);
                        os << x << " " << y << " " << psi[r][theta][phi].imag() << std::endl;
                    }
                }
            }
        }

    //    void verify_solution(const Parameters& params, const QuantumState& qs)
    //    {
    //        static constexpr complex_t CZERO{ 0.0, 0.0 };

    //        Array3d<complex_t> ratio(radial.size(), azimuthal.size(), polar.size());

    //        number_t alpha = params.h / (2.0 * params.mp);
    //        number_t E = params.ke / (qs.n * qs.n);

    //        for (std::size_t r = 1; r < radial.size() - 1; ++r)
    //        {
    //            for (std::size_t theta = 1; theta < polar.size() - 1; ++theta)
    //            {
    //                for (std::size_t phi = 1; phi < azimuthal.size() - 1; ++phi)
    //                {
    //                    if (psi[r][theta][phi] == CZERO)
    //                        continue;

    //                    number_t dV1 = radial[r] * radial[r] * radial.delta;
    //                    number_t dV2 = radial[r] * radial[r] * sin(polar[theta]) * sin(polar[theta]) * azimuthal.delta;
    //                    number_t dV3 = radial[r] * radial[r] * sin(polar[theta]) * polar.delta;

    //                    // number_t aS = r[i - 1] * r[i - 1] * alpha / (dV1 * r_p.delta);
    //                    // number_t aN = r[i] * r[i] * alpha / (dV1 * r_p.delta);
    //                    number_t r1 = radial[r] - radial.delta / 2.0;
    //                    number_t r2 = radial[r] + radial.delta / 2.0;

    //                    number_t aS = r1 * r1 * alpha / (dV1 * radial.delta);
    //                    number_t aN = r2 * r2 * alpha / (dV1 * radial.delta);
    //                    number_t aW = alpha / (dV2 * azimuthal.delta);
    //                    number_t aE = alpha / (dV2 * azimuthal.delta);
    //                    number_t aB = sin(polar[theta] - polar.delta / 2.0) * alpha / (dV3 * polar.delta);
    //                    number_t aT = sin(polar[theta] + polar.delta / 2.0) * alpha / (dV3 * polar.delta);

    //                    complex_t psi_s{ -aS * psi[r - 1][phi][theta] };
    //                    complex_t psi_n{ -aN * psi[r + 1][phi][theta] };
    //                    complex_t psi_w{ -aW * psi[r][phi - 1][theta] };
    //                    complex_t psi_e{ -aE * psi[r][phi + 1][theta] };
    //                    complex_t psi_b{ -aB * psi[r][phi][theta - 1] };
    //                    complex_t psi_t{ -aT * psi[r][phi][theta + 1] };

    //                    complex_t psi_p_k_rp{ -params.kp / radial[r] * psi[r][phi][theta] };

    //                    complex_t Fp{ psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + psi_p_k_rp };
    //                    number_t ap = aS + aN + aW + aE + aB + aT - E;
    //                    complex_t ap_psi_p{ -ap * psi[r][theta][phi] };

    //                    ratio[r][theta][phi] = Fp / ap_psi_p;
    //                }
    //            }
    //        }

    //        /* Compute error */
    //        number_t error = 0.0;
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
