#ifndef MANGH_HYDROGEN_SOLUTION_DATA_HPP
#define MANGH_HYDROGEN_SOLUTION_DATA_HPP

#include <fstream>

#include "numbers.hpp"
#include "linspace.hpp"
#include "array3d.hpp"
#include "Parameters.hpp"
#include "QuantumState.hpp"

namespace mangh::hydrogen
{
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

        /// @brief Azimuthal coordinate θ (angle in the XY-plane from the X-axis such that: 0 ≤ θ < 2π).
        linspace<number_t> theta;

        /// @brief Polar coordinate φ (angle from the positive Z axis such that: 0 ≤ φ ≤ π).
        linspace<number_t> phi;

        /// @brief Wave function (Ψ).
        array3d<complex_t> psi;

        /// @brief Maximum probability density value.
        number_t max_pd;

        /// @brief Maximum value of the real part of the wave function (Ψ).
        number_t max_real;

        //number_t E;
        number_t error;

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
         * @param radial - radial coordinate vector,
         * @param theta - azimuthal coordinate vector,
         * @param phis - polar coordinate vector.
         */
        Mesh(const linspace<number_t>& radial,
             const linspace<number_t>& theta,
             const linspace<number_t>& phi
            )
            : radial(radial)
            , theta(theta)
            , phi(phi)

            , psi(radial.size(), theta.size(), phi.size())

            , max_pd(std::numeric_limits<number_t>::min())
            , max_real(std::numeric_limits<number_t>::min())
            , error(std::numeric_limits<number_t>::min())
            , total_probability(std::numeric_limits<number_t>::min())
        {
            for(auto& th : this->theta) th += theta.delta / 2.0;
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        void compute_wave_function(const Parameters& param, const QuantumState& qs)
        {
            max_pd = std::numeric_limits<number_t>::min();
            max_real = std::numeric_limits<number_t>::min();

            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                // normalized radial distance from the nucleus:
                //      ρ = 2 * r / (n * a0):
                number_t rho = radial[i] / param.a0;

                for (std::size_t j = 0; j < theta.size(); ++j)
                {
                    for (std::size_t k = 0; k < phi.size(); ++k)
                    {
                        psi[i][j][k] = qs.WaveFunction(rho, theta[j], phi[k]);
                    
                        if (psi[i][j][k].real() > max_real) max_real = psi[i][j][k].real();
                    
                        number_t pd = norm(psi[i][j][k]);

                        if (pd > max_pd) max_pd = pd;
                    }
                }
            }
        }

        void normalize_psi()
        {
            static constexpr number_t M_PI = mangh::math::PI_V<number_t>;
            // static constexpr number_t M_PI = mangh::math::PI;

            total_probability = 0.0;

            /* Top pole */
            for (std::size_t i = 1; i < radial.size(); ++i)
            {
                number_t rlo = radial[i] - radial.delta / 2.0;
                number_t rhi = radial[i] + radial.delta / 2.0;
                number_t probability_per_vol =
                    2.0 / 3.0 * M_PI * (1 - cos(phi.delta / 2.0)) * (rhi * rhi * rhi - rlo * rlo * rlo) * norm(psi[i][0][0]);

                total_probability += probability_per_vol;
            }

            /* Central nodes */
            for (std::size_t i = 1; i < radial.size(); ++i)
            {
                for (std::size_t j = 0; j < theta.size(); ++j)
                {
                    for (std::size_t k = 1; k < phi.size() - 1; ++k)
                    {
                        number_t probability_per_vol = 
                            radial[i] * radial[i] * sin(phi[k]) * phi.delta * theta.delta * radial.delta * norm(psi[i][j][k]);
                        total_probability += probability_per_vol;
                    }
                }
            }

            /* Bottom pole */
            for (std::size_t i = 1; i < radial.size(); ++i)
            {
                number_t rlo = radial[i] - radial.delta / 2.0;
                number_t rhi = radial[i] + radial.delta / 2.0;
                number_t probability_per_vol = 
                    2.0 / 3.0 * M_PI * (1 - cos(phi.delta / 2.0)) * (rhi * rhi * rhi - rlo * rlo * rlo) * norm(psi[i][0][phi.size() - 1]);
                total_probability += probability_per_vol;
            }

            // auto integral_sqrt = sqrt(integral);
            // // auto integral_sqrt = integral / integral;

            // /* Normalize */
            // for (std::size_t i = 1; i < dd.n_r; ++i)
            // {
            //     for (std::size_t j = 0; j < dd.n_theta; ++j)
            //     {
            //         for (std::size_t k = 0; k < dd.n_phi; ++k)
            //         {
            //             psi[i][j][k] /= integral_sqrt;  // sqrt(integral);
            //         }
            //     }
            // }
        }

        void export_grid_data(std::string file_name/* , double max_pd, double max_real */)
        {
            std::ofstream myfile;
            myfile.open(file_name);

            myfile
                << radial.size() << " "
                << theta.size() << " "
                << phi.size() << " "
                << max_pd << " "
                << max_real
                << std::endl;

            myfile.close();
        }

        void export_probability_density(std::string file_name)
        {
            std::ofstream myfile_pd;
            myfile_pd.open(file_name);

            /* Export psi central nodes */
            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                for (std::size_t j = 0; j < 1; ++j)
                {
                    for (std::size_t k = 0; k < phi.size(); ++k)
                    {
                        number_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        number_t z = radial[i] * cos(phi[k]);
                        myfile_pd << x << " " << z << " " << norm(psi[i][j][k]) << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                for (std::size_t j = theta.size()/2; j < theta.size()/2 + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi.size(); ++k)
                    {
                        number_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        number_t z = radial[i] * cos(phi[k]);
                        myfile_pd << x << " " << z << " " << norm(psi[i][j][k]) << std::endl;
                    }
                }
            }

            myfile_pd.close();
        }

        void export_psi_real(std::string file_name)
        {
            std::ofstream myfile_real;
            myfile_real.open(file_name);

            std::size_t theta_slice = 0; // Should be a value between 0 and theta.size() / 2 - 2

            /* Export psi central nodes */
            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                for (std::size_t j = theta_slice; j < theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi.size(); ++k)
                    {
                        number_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        number_t z = radial[i] * cos(phi[k]);
                        myfile_real << x << " " << z << " " << psi[i][j][k].real() << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                for (std::size_t j = theta.size() / 2 + theta_slice; j < theta.size() / 2 + theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi.size(); ++k)
                    {
                        double x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        double z = radial[i] * cos(phi[k]);
                        myfile_real << x << " " << z << " " << psi[i][j][k].real() << std::endl;
                    }
                }
            }

            myfile_real.close();
        }

        void export_psi_im(std::string file_name)
        {
            std::ofstream myfile_im;
            myfile_im.open(file_name);

            std::size_t theta_slice = 0; // Should be a value between 0 and n_theta / 2 - 2

            /* Export psi central nodes */
            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                for (std::size_t j = theta_slice; j < theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi.size(); ++k)
                    {
                        number_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        number_t z = radial[i] * cos(phi[k]);
                        myfile_im << x << " " << z << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                for (std::size_t j = theta.size() / 2 + theta_slice; j < theta.size() / 2 + theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi.size(); ++k)
                    {
                        number_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        number_t z = radial[i] * cos(phi[k]);
                        myfile_im << x << " " << z << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            myfile_im.close();
        }

        void export_data_phi_2(std::string file_name)
        {
            std::ofstream myfile_phi_2;
            myfile_phi_2.open(file_name);

            /* Export psi central nodes */
            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                for (std::size_t j = 0; j < theta.size(); ++j)
                {
                    for (std::size_t k = phi.size() / 2; k < phi.size() / 2 + 1; ++k)
                    {
                        number_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        number_t y = radial[i] * sin(phi[k]) * sin(theta[j]);
                        myfile_phi_2 << x << " " << y << " " << psi[i][j][k].real() << std::endl;
                    }
                }
            }

            myfile_phi_2.close();
        }

        void export_data_phi_2_im(std::string file_name)
        {
            std::ofstream myfile_phi_2_im;
            myfile_phi_2_im.open(file_name);

            /* Export psi central nodes */
            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                for (std::size_t j = 0; j < theta.size(); ++j)
                {
                    for (std::size_t k = phi.size() / 2; k < phi.size() / 2 + 1; ++k)
                    {
                        number_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        number_t y = radial[i] * sin(phi[k]) * sin(theta[j]);
                        myfile_phi_2_im << x << " " << y << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            myfile_phi_2_im.close();
        }

        void verify_solution(const Parameters& params, const QuantumState& qs)
        {
            static constexpr complex_t CZERO{ 0.0, 0.0 };

            array3d<complex_t> ratio(radial.size(), theta.size(), phi.size());

            number_t alpha = params.h / (2.0 * params.mp);
            number_t E = params.ke / (qs.n * qs.n);

            for (std::size_t i = 1; i < radial.size() - 1; ++i)
            {
                for (std::size_t j = 1; j < theta.size() - 1; ++j)
                {
                    for (std::size_t k = 1; k < phi.size() - 1; ++k)
                    {
                        if (psi[i][j][k] == CZERO)
                            continue;

                        number_t dV1 = radial[i] * radial[i] * radial.delta;
                        number_t dV2 = radial[i] * radial[i] * sin(phi[k]) * sin(phi[k]) * theta.delta;
                        number_t dV3 = radial[i] * radial[i] * sin(phi[k]) * phi.delta;

                        // number_t aS = r[i - 1] * r[i - 1] * alpha / (dV1 * r_p.delta);
                        // number_t aN = r[i] * r[i] * alpha / (dV1 * r_p.delta);
                        number_t r1 = radial[i] - radial.delta / 2.0;
                        number_t r2 = radial[i] + radial.delta / 2.0;

                        number_t aS = r1 * r1 * alpha / (dV1 * radial.delta);
                        number_t aN = r2 * r2 * alpha / (dV1 * radial.delta);
                        number_t aW = alpha / (dV2 * theta.delta);
                        number_t aE = alpha / (dV2 * theta.delta);
                        number_t aB = sin(phi[k] - phi.delta / 2.0) * alpha / (dV3 * phi.delta);
                        number_t aT = sin(phi[k] + phi.delta / 2.0) * alpha / (dV3 * phi.delta);

                        complex_t psi_s{ -aS * psi[i - 1][j][k] };
                        complex_t psi_n{ -aN * psi[i + 1][j][k] };
                        complex_t psi_w{ -aW * psi[i][j - 1][k] };
                        complex_t psi_e{ -aE * psi[i][j + 1][k] };
                        complex_t psi_b{ -aB * psi[i][j][k - 1] };
                        complex_t psi_t{ -aT * psi[i][j][k + 1] };

                        complex_t psi_p_k_rp{ -params.kp / radial[i] * psi[i][j][k] };

                        complex_t Fp{ psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + psi_p_k_rp };
                        number_t ap = aS + aN + aW + aE + aB + aT - E;
                        complex_t ap_psi_p{ -ap * psi[i][j][k] };

                        ratio[i][j][k] = Fp / ap_psi_p;
                    }
                }
            }

            /* Compute error */
            error = 0.0;
            for (std::size_t i = 1; i < radial.size() - 1; ++i)
            {
                for (std::size_t j = 1; j < theta.size() - 1; ++j)
                {
                    for (std::size_t k = 1; k < phi.size() - 1; ++k)
                    {
                        if (fabs(ratio[i][j][k].real()) > 0.2)
                        {
                            error += fabs(1.0 - ratio[i][j][k].real());
                        }
                    }
                }
            }
            error = 100.0 * error / ((radial.size() - 2) * (theta.size() - 2) * (phi.size() - 2));
        }
    };
}

#endif /* !MANGH_HYDROGEN_SOLUTION_DATA_HPP */
