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

        /// @brief Radial coordinate (r primary).
        linspace<numeral_t> r_p;

        /// @brief Radial coordinate (r).
        linspace<numeral_t> r;

        /// @brief Polar angle (θ) (?)
        linspace<numeral_t> theta_p;

        /// @brief Azimuthal angle (ϕ primary)
        linspace<numeral_t> phi_p;

        /// @brief Azimuthal angle (ϕ)
        linspace<numeral_t> phi;

        /// @brief Wave function (Ψ).
        Array3d<complex_t> psi;

        /// @brief brief Probability density |Ψ|²
        Array3d<numeral_t> psi_square;

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
         * @param d - domain data.
         */
        SolutionData(
                const linspace<numeral_t>& radii, 
                const linspace<numeral_t>& thetas, 
                const linspace<numeral_t>& phis
            )
            : r_p(radii)
            , r(radii)
            , theta_p(thetas)
            , phi_p(phis)
            , phi(phis)     //, phi(dd.n_phi - 1)  // NOTE: n_phi - 1 !!!!!!!!!!!!!!!! Why?
            , psi(radii.size(), thetas.size(), phis.size())
            , psi_square(radii.size(), thetas.size(), phis.size())

            , max_pd(std::numeric_limits<numeral_t>::min())
            , max_real(std::numeric_limits<numeral_t>::min())
            , error(std::numeric_limits<numeral_t>::min())
            , integral(std::numeric_limits<numeral_t>::min())
        {
            // r & r_p vector 
            for(auto& radius : r) radius += radii.delta / 2.0;
            
            // theta_p vector
            for(auto& theta : theta_p) theta += thetas.delta / 2.0;

            //// phi vector
            //for (std::size_t i = 0; i < dd.n_phi - 1; ++i) { phi[i] = i * dphi + 0.5 * dphi; }

            //// phi_p vector
            //for (std::size_t i = 0; i < dd.n_phi; ++i) { phi_p[i] = i * dphi; }

            // phi & phi_p vector
            // for (std::size_t i = 0; i < phis.size(); ++i)
            // {
            //     phi_p[i] = i * dphi;
            //     phi[i] = i * dphi + 0.5 * dphi;
            // }
            for (auto& fi : phi) fi += phis.delta / 2.0;
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        void compute_wave_function(const PhysicalParameters& param, const QuantumState& qs)
        {
            max_pd = std::numeric_limits<numeral_t>::min();
            max_real = std::numeric_limits<numeral_t>::min();

            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                // normalized radial distance from the nucleus:
                //      ρ = 2 * r / (n * a0):
                numeral_t rho = r_p[i] / param.a0;

                for (std::size_t j = 0; j < theta_p.size(); ++j)
                {
                    for (std::size_t k = 0; k < phi_p.size(); ++k)
                    {
                        psi[i][j][k] = qs.WaveFunction(rho, theta_p[j], phi_p[k]);

                        psi_square[i][j][k] = norm(psi[i][j][k]);

                        if (psi_square[i][j][k] > max_pd)
                        {
                            max_pd = psi_square[i][j][k];
                        }

                        if (psi[i][j][k].real() > max_real)
                        {
                            max_real = psi[i][j][k].real();
                        }
                    }
                }
            }
        }

        void normalize_psi()
        {
            integral = 0.0;

            /* Top pole */
            for (std::size_t i = 1; i < r.size(); ++i)
            {
                //numeral_t R1 = r_p[i-1] + r_p.delta / 2.0;
                //numeral_t R2 = r_p[i] + r_p.delta / 2.0;
                numeral_t dV = 2.0 / 3.0 * M_PI * (1 - cos(phi.delta / 2.0)) * (r[i] * r[i] * r[i] - r[i-1] * r[i-1] * r[i-1]);
                numeral_t psi_sq = norm(psi[i][0][0]) * dV;
                integral += psi_sq;
            }

            /* Central nodes */
            for (std::size_t i = 1; i < r_p.size(); ++i)
            {
                for (std::size_t j = 0; j < theta_p.size(); ++j)
                {
                    for (std::size_t k = 1; k < phi_p.size() - 1; ++k)
                    {
                        numeral_t dV = r_p[i] * r_p[i] * sin(phi_p[k]) * phi_p.delta * theta_p.delta * r_p.delta;
                        numeral_t psi_sq = norm(psi[i][j][k]) * dV;
                        integral += psi_sq;
                    }
                }
            }

            /* Bottom pole */
            for (std::size_t i = 1; i < r_p.size(); ++i)
            {
                numeral_t dV = 2.0 / 3.0 * M_PI * (1 - cos(phi_p.delta / 2.0)) * (r[i] * r[i] * r[i] - r[i-1] * r[i-1] * r[i-1]);
                numeral_t psi_sq = norm(psi[i][0][phi_p.size() - 1]) * dV;
                integral += psi_sq;
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

        // void compute_probability_density()
        // {
        //     for (std::size_t i = 0; i < dd.n_r; ++i)
        //     {
        //         for (std::size_t j = 0; j < dd.n_theta; ++j)
        //         {
        //             for (std::size_t k = 0; k < dd.n_phi; ++k)
        //             {
        //                 psi_square[i][j][k] = norm(psi[i][j][k]);
        //             }
        //         }
        //     }
        // }

        // numeral_t compute_max_pd()
        // {
        //     //numeral_t max_pd = -1e+8;
        //     numeral_t max_pd = std::numeric_limits<numeral_t>::min();
        //     for (std::size_t i = 0; i < dd.n_r; ++i)
        //     {
        //         for (std::size_t j = 0; j < dd.n_theta; ++j)
        //         {
        //             for (std::size_t k = 0; k < dd.n_phi; ++k)
        //             {
        //                 if (psi_square[i][j][k] > max_pd)
        //                 {
        //                     max_pd = psi_square[i][j][k];
        //                 }
        //             }
        //         }
        //     }
        //     return max_pd;
        // }

        // numeral_t compute_max_psi_real()
        // {
        //     //numeral_t max_real = -1e+8;
        //     numeral_t max_real = std::numeric_limits<numeral_t>::min();
        //     for (std::size_t i = 0; i < dd.n_r; ++i)
        //     {
        //         for (std::size_t j = 0; j < dd.n_theta; ++j)
        //         {
        //             for (std::size_t k = 0; k < dd.n_phi; ++k)
        //             {
        //                 if (psi[i][j][k].real() > max_real)
        //                 {
        //                     max_real = psi[i][j][k].real();
        //                 }
        //             }
        //         }
        //     }
        //     return max_real;
        // }

        //////////////////////////////////////////////////////////////////////////////////////////////
        // void export_psi(
        //     const std::string& file_name,
        //     const PhysicalParameters& param,
        //     const QuantumState& qs,
        //     const linspace<numeral_t>& radii, 
        //     const linspace<numeral_t>& thetas, 
        //     const linspace<numeral_t>& phis
        //     )
        // {

        //     std::ofstream myfile;
        //     myfile.open(file_name);

        //     constexpr auto max_precision{std::numeric_limits<numeral_t>::digits10 + 1};     
        //     // myfile.imbue(std::locale(myfile.getloc(), new punct_facet<char, ','>));

        //     // myfile
        //     //     << std::setprecision(max_precision)
        //     //     << "radius[0] = " << radii[0] << "\t"
        //     //     << "radius[delta] = " << radii.delta << "\t"
        //     //     << "radius[last] = " << radii[radii.size() - 1]
        //     //     << std::endl
        //     //     << "thetas[0] = " << thetas[0] << "\t"
        //     //     << "thetas[delta] = " << thetas.delta << "\t"
        //     //     << "thetas[last] = " << thetas[thetas.size() - 1]
        //     //     << std::endl
        //     //     << "phis[0] = " << phis[0] << "\t"
        //     //     << "phis[delta] = " << phis.delta << "\t"
        //     //     << "phis[last] = " << phis[phis.size() - 1]
        //     //     << std::endl
        //     //     << std::endl;

        //     myfile
        //         // << "r\t"
        //         // << "theta\t"
        //         // << "phi\t"
        //         << "x\t"
        //         << "y\t"
        //         << "z\t"
        //         // << "psi_re\t"
        //         // << "psi_im\t"
        //         << "psi2"
        //         // << std::scientific
        //         << std::setprecision(max_precision)
        //         << std::endl;

        //     // for (const auto& radius : radii)
        //     for (size_t r = 0; r < radii.size(); ++r)
        //     {
        //         numeral_t radius = radii[r];
        //         // normalized radial distance from the nucleus: ρ = 2 * radius / (n * a0):
        //         // numeral_t rho = radius / (qs.n * param.a0);
        //         numeral_t rho = radius / param.a0;

        //         // for (const auto& theta : thetas)
        //         for (size_t t = 0; t < thetas.size(); ++t)
        //         {
        //             numeral_t theta = thetas[t];
        //             // for (const auto& phi : phis)
        //             for (size_t p = 0; p < phis.size(); ++p)
        //             {
        //                 numeral_t phi = phis[p];
        //                 complex_t PSI = psi[r][t][p] = qs.WaveFunction(rho, theta, phi);

        //                 // https://en.wikipedia.org/wiki/Spherical_coordinate_system

        //                 // (mathematics convention)
        //                 numeral_t x = radius * std::sin(phi) * std::cos(theta);
        //                 numeral_t y = radius * std::sin(phi) * std::sin(theta);
        //                 numeral_t z = radius * std::cos(phi);

        //                 // // (physics convention)
        //                 // numeral_t x = radius * std::sin(theta) * std::cos(phi);
        //                 // numeral_t y = radius * std::sin(theta) * std::sin(phi);
        //                 // numeral_t z = radius * std::cos(theta)

        //                 myfile
        //                     // << radius << "\t"
        //                     // << theta << "\t"
        //                     // << phi << "\t"
        //                     << x << "\t"
        //                     << y << "\t"
        //                     << z << "\t"
        //                     // << PSI << "\t"
        //                     // << PSI.real() << "\t"
        //                     // << PSI.imag() << "\t"
        //                     << norm(PSI)
        //                     << std::endl;
        //             }
        //         }
        //     }
        //     myfile.close();
        // }

        void export_grid_data(std::string file_name/* , double max_pd, double max_real */)
        {
            std::ofstream myfile;
            myfile.open(file_name);

            myfile
                << r_p.size() << " "
                << theta_p.size() << " "
                << phi_p.size() << " "
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
            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                for (std::size_t j = 0; j < 1; ++j)
                {
                    for (std::size_t k = 0; k < phi_p.size(); ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t z = r_p[i] * cos(phi_p[k]);
                        myfile_pd << x << " " << z << " " << psi_square[i][j][k] << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                for (std::size_t j = theta_p.size()/2; j < theta_p.size()/2 + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi_p.size(); ++k)
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
            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                for (std::size_t j = theta_slice; j < theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi_p.size(); ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t z = r_p[i] * cos(phi_p[k]);
                        myfile_real << x << " " << z << " " << psi[i][j][k].real() << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                for (std::size_t j = theta_p.size() / 2 + theta_slice; j < theta_p.size() / 2 + theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi_p.size(); ++k)
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
            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                for (std::size_t j = theta_slice; j < theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi_p.size(); ++k)
                    {
                        numeral_t x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        numeral_t z = r_p[i] * cos(phi_p[k]);
                        myfile_im << x << " " << z << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                for (std::size_t j = theta_p.size() / 2 + theta_slice; j < theta_p.size() / 2 + theta_slice + 1; ++j)
                {
                    for (std::size_t k = 0; k < phi_p.size(); ++k)
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
            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                for (std::size_t j = 0; j < theta_p.size(); ++j)
                {
                    for (std::size_t k = phi_p.size() / 2; k < phi_p.size() / 2 + 1; ++k)
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
            for (std::size_t i = 0; i < r_p.size(); ++i)
            {
                for (std::size_t j = 0; j < theta_p.size(); ++j)
                {
                    for (std::size_t k = phi_p.size() / 2; k < phi_p.size() / 2 + 1; ++k)
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
            Array3d<complex_t> ratio(r_p.size(), theta_p.size(), phi_p.size());

            numeral_t alpha = params.h / (2.0 * params.mp);
            numeral_t E = params.ke / (qs.n * qs.n);

            for (std::size_t i = 1; i < r_p.size() - 1; ++i)
            {
                for (std::size_t j = 1; j < theta_p.size() - 1; ++j)
                {
                    for (std::size_t k = 1; k < phi_p.size() - 1; ++k)
                    {
                        if (psi[i][j][k] == CZERO)
                            continue;

                        numeral_t dV1 = r_p[i] * r_p[i] * r_p.delta;
                        numeral_t dV2 = r_p[i] * r_p[i] * sin(phi_p[k]) * sin(phi_p[k]) * theta_p.delta;
                        numeral_t dV3 = r_p[i] * r_p[i] * sin(phi_p[k]) * phi_p.delta;

                        numeral_t aS = r[i - 1] * r[i - 1] * alpha / (dV1 * r_p.delta);
                        numeral_t aN = r[i] * r[i] * alpha / (dV1 * r_p.delta);
                        numeral_t aW = alpha / (dV2 * theta_p.delta);
                        numeral_t aE = alpha / (dV2 * theta_p.delta);
                        numeral_t aB = sin(phi[k - 1]) * alpha / (dV3 * phi_p.delta);
                        numeral_t aT = sin(phi[k]) * alpha / (dV3 * phi_p.delta);

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
            for (std::size_t i = 1; i < r_p.size() - 1; ++i)
            {
                for (std::size_t j = 1; j < theta_p.size() - 1; ++j)
                {
                    for (std::size_t k = 1; k < phi_p.size() - 1; ++k)
                    {
                        if (fabs(ratio[i][j][k].real()) > 0.2)
                        {
                            error += fabs(1.0 - ratio[i][j][k].real());
                        }
                    }
                }
            }
            error = 100.0 * error / ((r_p.size() - 2) * (theta_p.size() - 2) * (phi_p.size() - 2));
        }
    };
}

#endif /* !SOLUTION_DATA_HPP */
