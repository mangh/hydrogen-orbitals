#ifndef SOLUTION_DATA_HPP
#define SOLUTION_DATA_HPP

#include <fstream>

#include "Array3d.hpp"
#include "DomainData.hpp"
#include "PhysicalParameters.hpp"
#include "QuantumState.hpp"

namespace HydrogenOrbitals
{
    constexpr numeral_t M_PI = 3.141592653589793238462643383279502884L;
    constexpr complex_t CZERO{ 0.0, 0.0 };

    /**
     * @brief Solution data.
     *
     */
    struct SolutionData
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constants
        ///

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief brief description
        const DomainData& dd;

        /// @brief Radial coordinate (r).
        std::vector<numeral_t> r;

        /// @brief brief description
        std::vector<numeral_t> r_p;

        /// @brief Polar angle (θ) (?)
        std::vector<numeral_t> theta_p;

        /// @brief Azimuthal angle (ϕ) (?)
        std::vector<numeral_t> phi;

        /// @brief brief description
        std::vector<numeral_t> phi_p;

        /// @brief Wave function (Ψ).
        Array3d<complex_t> psi;

        /// @brief brief Probability density |Ψ|²
        Array3d<numeral_t> psi_square;

        /// @brief brief description
        numeral_t E;

        /// @brief brief description
        numeral_t error;

        // Normalizing factor
        numeral_t integral;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        SolutionData() = delete;

        /**
         * @brief SolutionData constructor.
         *
         * @param d - domain data.
         */
        SolutionData(const DomainData& domain_data)
            : dd(domain_data)
            , r(dd.n_r)
            , r_p(dd.n_r)
            , theta_p(dd.n_theta)
            , phi(dd.n_phi - 1)  // NOTE: n_phi - 1 !!!!!!!!!!!!!!!! Why?
            , phi_p(dd.n_phi)
            , psi(dd)
            , psi_square(dd)
        {
            numeral_t dr = dd.R / dd.n_r;
            numeral_t dtheta = 2.0 * M_PI / dd.n_theta;
            numeral_t dphi = M_PI / (dd.n_phi - 1);

            // r & r_p vector 
            for (std::size_t i = 0; i < dd.n_r; ++i) { r_p[i] = i * dr; r[i] = r_p[i] + 0.5 * dr; }

            // theta_p vector
            for (std::size_t i = 0; i < dd.n_theta; ++i) { theta_p[i] = i * dtheta + 0.5 * dtheta; }

            // phi vector
            for (std::size_t i = 0; i < dd.n_phi - 1; ++i) { phi[i] = i * dphi + 0.5 * dphi; }

            // phi_p vector
            for (std::size_t i = 0; i < dd.n_phi; ++i) { phi_p[i] = i * dphi; }
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        void compute_wave_function(const PhysicalParameters& param, const QuantumState& qs)
        {
            for (std::size_t i = 0; i < 1; ++i)
            {
                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        psi[i][j][k] = CZERO;
                    }
                }
            }

            for (std::size_t i = 1; i < dd.n_r; ++i)
            {
                // normalized radial distance from the nucleus: ρ = 2 * r / (n * a0):
                // numeral_t rho = r_p[i] / (qs.n * param.a0);
                numeral_t rho = r_p[i] / param.a0;

                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        psi[i][j][k] = qs.WaveFunction(rho, theta_p[j], phi_p[k]);
                    }
                }
            }
        }

        void normalize_psi()
        {
            // numeral_t integral{ 0.0 };
            integral = 0.0;

            numeral_t dr = dd.R / dd.n_r;
            numeral_t dtheta = 2.0 * M_PI / dd.n_theta;
            numeral_t dphi = M_PI / (dd.n_phi - 1);

            /* Top pole */
            for (std::size_t i = 1; i < dd.n_r; ++i)
            {
                numeral_t dV = 2.0 / 3.0 * M_PI * (1 - cos(dphi / 2)) * (r[i] * r[i] * r[i] - r[i - 1] * r[i - 1] * r[i - 1]);
                numeral_t psi_sq = norm(psi[i][0][0]) * dV;
                integral += psi_sq;
            }

            /* Central nodes */
            for (std::size_t i = 1; i < dd.n_r; ++i)
            {
                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = 1; k < dd.n_phi - 1; ++k)
                    {
                        numeral_t dV = r_p[i] * r_p[i] * sin(phi_p[k]) * dphi * dtheta * dr;
                        numeral_t psi_sq = norm(psi[i][j][k]) * dV;
                        integral += psi_sq;
                    }
                }
            }

            /* Bottom pole */
            for (std::size_t i = 1; i < dd.n_r; ++i)
            {
                numeral_t dV = 2.0 / 3.0 * M_PI * (1 - cos(dphi / 2.0)) * (r[i] * r[i] * r[i] - r[i - 1] * r[i - 1] * r[i - 1]);
                numeral_t psi_sq = norm(psi[i][0][dd.n_phi - 1]) * dV;
                integral += psi_sq;
            }

            auto integral_sqrt = sqrt(integral);
            // auto integral_sqrt = integral / integral;

            /* Normalize */
            for (std::size_t i = 1; i < dd.n_r; ++i)
            {
                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        psi[i][j][k] /= integral_sqrt;  // sqrt(integral);
                    }
                }
            }
        }

        void compute_probability_density()
        {
            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        psi_square[i][j][k] = norm(psi[i][j][k]);
                    }
                }
            }
        }

        numeral_t compute_max_pd()
        {
            numeral_t max_pd = -1e+8;
            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        if (psi_square[i][j][k] > max_pd)
                        {
                            max_pd = psi_square[i][j][k];
                        }
                    }
                }
            }
            return max_pd;
        }

        numeral_t compute_max_psi_real()
        {
            numeral_t max_real = -1e+8;
            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        if (psi[i][j][k].real() > max_real)
                        {
                            max_real = psi[i][j][k].real();
                        }
                    }
                }
            }
            return max_real;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////
        void export_grid_data(std::string file_name, double max_pd, double max_real)
        {
            std::ofstream myfile;
            myfile.open(file_name);

            myfile
                << dd.n_r << " "
                << dd.n_theta << " "
                << dd.n_phi << " "
                << max_pd << " "
                << max_real << std::endl;

            myfile.close();
        }

        void export_probability_density(std::string file_name)
        {
            std::ofstream myfile_pd;
            myfile_pd.open(file_name);

            /* Export psi central nodes */
            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = 0; j < 1; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t z = r_p[i] * cos(phi_p[k]);
                        myfile_pd << x << " " << z << " " << psi_square[i][j][k] << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = dd.n_theta/2; j < dd.n_theta/2 + 1; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t z = r_p[i] * cos(phi_p[k]);
                        myfile_pd << x << " " << z << " " << psi_square[i][j][k] << std::endl;
                    }
                }
            }

            myfile_pd.close();
        }

        void export_psi_real(std::string file_name)
        {
            std::ofstream myfile_real;
            myfile_real.open(file_name);

            std::size_t theta_slice = 0; // Should be a value between 0 and n_theta / 2 - 2

            /* Export psi central nodes */
            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = theta_slice; j < theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t z = r_p[i] * cos(phi_p[k]);
                        myfile_real << x << " " << z << " " << psi[i][j][k].real() << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = dd.n_theta / 2 + theta_slice; j < dd.n_theta / 2 + theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        double x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        double z = r_p[i] * cos(phi_p[k]);
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
            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = theta_slice; j < theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t z = r_p[i] * cos(phi_p[k]);
                        myfile_im << x << " " << z << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = dd.n_theta / 2 + theta_slice; j < dd.n_theta / 2 + theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < dd.n_phi; ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t z = r_p[i] * cos(phi_p[k]);
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
            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = dd.n_phi / 2; k < dd.n_phi / 2 + 1; ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t y = r_p[i] * sin(phi_p[k]) * sin(theta_p[j]);
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
            for (std::size_t i = 0; i < dd.n_r; ++i)
            {
                for (std::size_t j = 0; j < dd.n_theta; ++j)
                {
                    for (std::size_t k = dd.n_phi / 2; k < dd.n_phi / 2 + 1; ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t y = r_p[i] * sin(phi_p[k]) * sin(theta_p[j]);
                        myfile_phi_2_im << x << " " << y << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            myfile_phi_2_im.close();
        }

        void verify_solution(const PhysicalParameters& params, const QuantumState& qs)
        {
            Array3d<complex_t> ratio(dd);

            numeral_t dr = dd.R / dd.n_r;
            numeral_t dtheta = 2.0 * M_PI / dd.n_theta;
            numeral_t dphi = M_PI / (dd.n_phi - 1);
            numeral_t alpha = params.h / (2.0 * params.mp);
            numeral_t E = params.ke / (qs.n * qs.n);

            for (std::size_t i = 1; i < dd.n_r - 1; ++i)
            {
                for (std::size_t j = 1; j < dd.n_theta - 1; ++j)
                {
                    for (std::size_t k = 1; k < dd.n_phi - 1; ++k)
                    {
                        if (psi[i][j][k] == CZERO)
                            continue;

                        numeral_t dV1 = r_p[i] * r_p[i] * dr;
                        numeral_t dV2 = r_p[i] * r_p[i] * sin(phi_p[k]) * sin(phi_p[k]) * dtheta;
                        numeral_t dV3 = r_p[i] * r_p[i] * sin(phi_p[k]) * dphi;

                        numeral_t aS = r[i - 1] * r[i - 1] * alpha / (dV1 * dr);
                        numeral_t aN = r[i] * r[i] * alpha / (dV1 * dr);
                        numeral_t aW = alpha / (dV2 * dtheta);
                        numeral_t aE = alpha / (dV2 * dtheta);
                        numeral_t aB = sin(phi[k - 1]) * alpha / (dV3 * dphi);
                        numeral_t aT = sin(phi[k]) * alpha / (dV3 * dphi);

                        complex_t psi_s{ -aS * psi[i - 1][j][k] };
                        complex_t psi_n{ -aN * psi[i + 1][j][k] };
                        complex_t psi_w{ -aW * psi[i][j - 1][k] };
                        complex_t psi_e{ -aE * psi[i][j + 1][k] };
                        complex_t psi_b{ -aB * psi[i][j][k - 1] };
                        complex_t psi_t{ -aT * psi[i][j][k + 1] };

                        complex_t psi_p_k_rp{ -params.kp / r_p[i] * psi[i][j][k] };

                        complex_t Fp{ psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + psi_p_k_rp };
                        numeral_t ap = aS + aN + aW + aE + aB + aT - E;
                        complex_t ap_psi_p{ -ap * psi[i][j][k] };

                        ratio[i][j][k] = Fp / ap_psi_p;
                    }
                }
            }

            /* Compute error */
            error = 0.0;
            for (std::size_t i = 1; i < dd.n_r - 1; ++i)
            {
                for (std::size_t j = 1; j < dd.n_theta - 1; ++j)
                {
                    for (std::size_t k = 1; k < dd.n_phi - 1; ++k)
                    {
                        if (fabs(ratio[i][j][k].real()) > 0.2)
                        {
                            error += fabs(1.0 - ratio[i][j][k].real());
                        }
                    }
                }
            }
            error = 100.0 * error / ((dd.n_r - 2) * (dd.n_theta - 2) * (dd.n_phi - 2));
        }
    };
}

#endif /* !SOLUTION_DATA_HPP */
