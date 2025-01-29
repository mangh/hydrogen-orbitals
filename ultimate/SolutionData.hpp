#ifndef SOLUTION_DATA_HPP
#define SOLUTION_DATA_HPP

#include <fstream>

#include "linspace.hpp"
#include "Numeral.hpp"

#include "Array3d.hpp"
#include "PhysicalParameters.hpp"
#include "QuantumState.hpp"

template <class charT, charT sep>
class punct_facet: public std::numpunct<charT>
{
protected:
    charT do_decimal_point() const { return sep; }
};

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

        /// @brief Radial coordinate.
        linspace<numeral_t> radial;

        /// @brief Azimuthal coordinate θ (angle in the XY-plane from the X-axis such that: 0 ≤ θ < 2π).
        linspace<numeral_t> theta;

        /// @brief Polar coordinate φ (angle from the positive Z axis such that: 0 ≤ ϕ ≤ π).
        linspace<numeral_t> phi;

        /// @brief Wave function (Ψ).
        Array3d<complex_t> psi;

        numeral_t max_pd;
        numeral_t max_real;
        //numeral_t E;
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
         * @param radial - radial coordinate vector,
         * @param theta - azimuthal coordinate vector,
         * @param phis - polar coordinate vector.
         */
        SolutionData(
                const linspace<numeral_t>& radial, 
                const linspace<numeral_t>& theta, 
                const linspace<numeral_t>& phi
            )
            : radial(radial)
            , theta(theta)
            , phi(phi)

            , psi(radial.size(), theta.size(), phi.size())

            , max_pd(std::numeric_limits<numeral_t>::min())
            , max_real(std::numeric_limits<numeral_t>::min())
            , error(std::numeric_limits<numeral_t>::min())
            , integral(std::numeric_limits<numeral_t>::min())
        {
            for(auto& th : this->theta) th += theta.delta / 2.0;
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        void compute_wave_function(const PhysicalParameters& param, const QuantumState& qs)
        {
            max_pd = std::numeric_limits<numeral_t>::min();
            max_real = std::numeric_limits<numeral_t>::min();

            for (std::size_t i = 0; i < radial.size(); ++i)
            {
                // normalized radial distance from the nucleus:
                //      ρ = 2 * r / (n * a0):
                numeral_t rho = radial[i] / param.a0;

                for (std::size_t j = 0; j < theta.size(); ++j)
                {
                    for (std::size_t k = 0; k < phi.size(); ++k)
                    {
                        psi[i][j][k] = qs.WaveFunction(rho, theta[j], phi[k]);
                    
                        if (psi[i][j][k].real() > max_real) max_real = psi[i][j][k].real();
                    
                        numeral_t pd = norm(psi[i][j][k]);

                        if (pd > max_pd) max_pd = pd;
                    }
                }
            }
        }

        void normalize_psi()
        {
            integral = 0.0;

            /* Top pole */
            for (std::size_t i = 1; i < radial.size(); ++i)
            {
                // d = r.delta / 2.0
                // r = r[i]
                // r1 = r[i] - d
                // r2 = r[i] + d
                // dvol = r2^3 - r1^3 = (r + d)^3 - (r - d)^3
                //      = 2 * d * (3*r^2 + d^2)
                numeral_t dvol = radial.delta * (3.0 * radial[i] * radial[i] + radial.delta * radial.delta / 4.0);
                numeral_t probability_per_vol =
                    2.0 / 3.0 * M_PI * (1 - cos(phi.delta / 2.0)) * dvol * norm(psi[i][0][0]);

                integral += probability_per_vol;
            }

            /* Central nodes */
            for (std::size_t i = 1; i < radial.size(); ++i)
            {
                for (std::size_t j = 0; j < theta.size(); ++j)
                {
                    for (std::size_t k = 1; k < phi.size() - 1; ++k)
                    {
                        numeral_t probability_per_vol = 
                            radial[i] * radial[i] * sin(phi[k]) * phi.delta * theta.delta * radial.delta * norm(psi[i][j][k]);
                        integral += probability_per_vol;
                    }
                }
            }

            /* Bottom pole */
            for (std::size_t i = 1; i < radial.size(); ++i)
            {
                numeral_t dVol = radial.delta * (3.0 * radial[i] * radial[i] + radial.delta * radial.delta / 4.0);
                numeral_t probability_per_vol = 
                    2.0 / 3.0 * M_PI * (1 - cos(phi.delta / 2.0)) * dVol * norm(psi[i][0][phi.size() - 1]);
                integral += probability_per_vol;
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
                        numeral_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        numeral_t z = radial[i] * cos(phi[k]);
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
                        numeral_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        numeral_t z = radial[i] * cos(phi[k]);
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
                        numeral_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        numeral_t z = radial[i] * cos(phi[k]);
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
                        numeral_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        numeral_t z = radial[i] * cos(phi[k]);
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
                        numeral_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        numeral_t z = radial[i] * cos(phi[k]);
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
                        numeral_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        numeral_t y = radial[i] * sin(phi[k]) * sin(theta[j]);
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
                        numeral_t x = radial[i] * sin(phi[k]) * cos(theta[j]);
                        numeral_t y = radial[i] * sin(phi[k]) * sin(theta[j]);
                        myfile_phi_2_im << x << " " << y << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            myfile_phi_2_im.close();
        }

        void verify_solution(const PhysicalParameters& params, const QuantumState& qs)
        {
            Array3d<complex_t> ratio(radial.size(), theta.size(), phi.size());

            numeral_t alpha = params.h / (2.0 * params.mp);
            numeral_t E = params.ke / (qs.n * qs.n);

            for (std::size_t i = 1; i < radial.size() - 1; ++i)
            {
                for (std::size_t j = 1; j < theta.size() - 1; ++j)
                {
                    for (std::size_t k = 1; k < phi.size() - 1; ++k)
                    {
                        if (psi[i][j][k] == CZERO)
                            continue;

                        numeral_t dV1 = radial[i] * radial[i] * radial.delta;
                        numeral_t dV2 = radial[i] * radial[i] * sin(phi[k]) * sin(phi[k]) * theta.delta;
                        numeral_t dV3 = radial[i] * radial[i] * sin(phi[k]) * phi.delta;

                        // numeral_t aS = r[i - 1] * r[i - 1] * alpha / (dV1 * r_p.delta);
                        // numeral_t aN = r[i] * r[i] * alpha / (dV1 * r_p.delta);
                        numeral_t r1 = radial[i] - radial.delta / 2.0;
                        numeral_t r2 = radial[i] + radial.delta / 2.0;

                        numeral_t aS = r1 * r1 * alpha / (dV1 * radial.delta);
                        numeral_t aN = r2 * r2 * alpha / (dV1 * radial.delta);
                        numeral_t aW = alpha / (dV2 * theta.delta);
                        numeral_t aE = alpha / (dV2 * theta.delta);
                        numeral_t aB = sin(phi[k] - phi.delta / 2.0) * alpha / (dV3 * phi.delta);
                        numeral_t aT = sin(phi[k] + phi.delta / 2.0) * alpha / (dV3 * phi.delta);

                        complex_t psi_s{ -aS * psi[i - 1][j][k] };
                        complex_t psi_n{ -aN * psi[i + 1][j][k] };
                        complex_t psi_w{ -aW * psi[i][j - 1][k] };
                        complex_t psi_e{ -aE * psi[i][j + 1][k] };
                        complex_t psi_b{ -aB * psi[i][j][k - 1] };
                        complex_t psi_t{ -aT * psi[i][j][k + 1] };

                        complex_t psi_p_k_rp{ -params.kp / radial[i] * psi[i][j][k] };

                        complex_t Fp{ psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + psi_p_k_rp };
                        numeral_t ap = aS + aN + aW + aE + aB + aT - E;
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

#endif /* !SOLUTION_DATA_HPP */
