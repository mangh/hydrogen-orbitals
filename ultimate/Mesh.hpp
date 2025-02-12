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

        /// @brief Polar coordinate 0 ≤ θ ≤ π (angle from the positive Z axis).
        linspace<number_t> theta;

        /// @brief Azimuthal coordinate 0 ≤ ϕ < 2π (angle from the X-axis in the XY-plane).
        linspace<number_t> phi;

        /// @brief Wave function (Ψ).
        array3d<complex_t> psi;

        /// @brief Maximum probability density value.
        number_t max_pd;

        /// @brief Maximum value of the real part of the wave function (Ψ).
        number_t max_real;

        //number_t error;

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
         * @param theta - polar coordinate vector,
         * @param phi - azimuthal coordinate vector.
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
            //, error(std::numeric_limits<number_t>::min())
            , total_probability(std::numeric_limits<number_t>::min())
        {
            for(auto& p : this->phi) p += phi.delta / 2.0;
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        void compute_wave_function(const QuantumState& qs)
        {
            max_pd = std::numeric_limits<number_t>::min();
            max_real = std::numeric_limits<number_t>::min();
            complex_t c;
            number_t pd;

            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = 0; it < theta.size(); ++it)
                {
                    for (std::size_t ip = 0; ip < phi.size(); ++ip)
                    {
                        psi[ir][it][ip] = c = qs.WaveFunction(radial[ir], theta[it], phi[ip]);

                        if (c.real() > max_real)
                            max_real = c.real();

                        if ((pd = norm(c)) > max_pd)
                            max_pd = pd;
                    }
                }
            }
        }

        void normalize_psi()
        {
            static constexpr number_t M_PI = mangh::math::PI_V<number_t>;

            total_probability = 0.0;

            /* Top pole */
            for (std::size_t ir = 1; ir < radial.size(); ++ir)
            {
                number_t rlo = radial[ir] - radial.delta / 2.0;
                number_t rhi = radial[ir] + radial.delta / 2.0;
                number_t probability_per_vol =
                    2.0 / 3.0 * M_PI * (1 - cos(theta.delta / 2.0)) * (rhi * rhi * rhi - rlo * rlo * rlo) * norm(psi[ir][0][0]);

                total_probability += probability_per_vol;
            }

            /* Central nodes */
            for (std::size_t ir = 1; ir < radial.size(); ++ir)
            {
                for (std::size_t it = 1; it < theta.size(); ++it)
                {
                    for (std::size_t ip = 0; ip < phi.size() - 1; ++ip)
                    {
                        number_t probability_per_vol = 
                            radial[ir] * radial[ir] * sin(theta[it]) * theta.delta * phi.delta * radial.delta * norm(psi[ir][it][ip]);
                        total_probability += probability_per_vol;
                    }
                }
            }

            /* Bottom pole */
            for (std::size_t ir = 1; ir < radial.size(); ++ir)
            {
                number_t rlo = radial[ir] - radial.delta / 2.0;
                number_t rhi = radial[ir] + radial.delta / 2.0;
                number_t probability_per_vol = 
                    2.0 / 3.0 * M_PI * (1 - cos(theta.delta / 2.0)) * (rhi * rhi * rhi - rlo * rlo * rlo) * norm(psi[ir][theta.size() - 1][0]);
                total_probability += probability_per_vol;
            }

            ///* Normalize */
            //auto integral_sqrt = sqrt(total_probability);
            //for (std::size_t ir = 1; ir < radial.size(); ++ir)
            //{
            //    for (std::size_t ith = 0; ith < theta.size(); ++ith)
            //    {
            //        for (std::size_t iph = 0; iph < phi.size(); ++iph)
            //        {
            //            psi[ir][ith][iph] /= integral_sqrt;
            //        }
            //    }
            //}
        }

        void export_grid_data(std::string file_name) const
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
            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = 0; it < theta.size(); ++it)
                {
                    for (std::size_t ip = 0; ip < 1; ++ip)
                    {
                        number_t x = radial[ir] * sin(theta[it]) * cos(phi[ip]);
                        number_t z = radial[ir] * cos(theta[it]);
                        myfile_pd << x << " " << z << " " << norm(psi[ir][it][ip]) << std::endl;
                    }
                }
            }

            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = 0; it < theta.size(); ++it)
                {
                    for (std::size_t ip = phi.size() / 2; ip < phi.size() / 2 + 1; ++ip)
                    {
                        number_t x = radial[ir] * sin(theta[it]) * cos(phi[ip]);
                        number_t z = radial[ir] * cos(theta[it]);
                        myfile_pd << x << " " << z << " " << norm(psi[ir][it][ip]) << std::endl;
                    }
                }
            }

            myfile_pd.close();
        }

        void export_psi_real(std::string file_name)
        {
            std::ofstream myfile_real;
            myfile_real.open(file_name);

            std::size_t phi_slice = 0; // Should be a value between 0 and phi.size() / 2 - 2

            /* Export psi central nodes */
            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = 0; it < theta.size(); ++it)
                {
                    for (std::size_t ip = phi_slice; ip < phi_slice + 1; ++ip)
                    {
                        number_t x = radial[ir] * sin(theta[it]) * cos(phi[ip]);
                        number_t z = radial[ir] * cos(theta[it]);
                        myfile_real << x << " " << z << " " << psi[ir][it][ip].real() << std::endl;
                    }
                }
            }

            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = 0; it < theta.size(); ++it)
                {
                    for (std::size_t ip = phi.size() / 2 + phi_slice; ip < phi.size() / 2 + phi_slice + 1; ++ip)
                    {
                        double x = radial[ir] * sin(theta[it]) * cos(phi[ip]);
                        double z = radial[ir] * cos(theta[it]);
                        myfile_real << x << " " << z << " " << psi[ir][it][ip].real() << std::endl;
                    }
                }
            }

            myfile_real.close();
        }

        void export_psi_im(std::string file_name)
        {
            std::ofstream myfile_im;
            myfile_im.open(file_name);

            std::size_t phi_slice = 0; // Should be a value between 0 and theta.size() / 2 - 2

            /* Export psi central nodes */
            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = 0; it < theta.size(); ++it)
                {
                    for (std::size_t ip = phi_slice; ip < phi_slice + 1; ++ip)
                    {
                        number_t x = radial[ir] * sin(theta[it]) * cos(phi[ip]);
                        number_t z = radial[ir] * cos(theta[it]);
                        myfile_im << x << " " << z << " " << psi[ir][it][ip].imag() << std::endl;
                    }
                }
            }

            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = 0; it < theta.size(); ++it)
                {
                    for (std::size_t ip = phi.size() / 2 + phi_slice; ip < phi.size() / 2 + phi_slice + 1; ++ip)
                    {
                        number_t x = radial[ir] * sin(theta[it]) * cos(phi[ip]);
                        number_t z = radial[ir] * cos(theta[it]);
                        myfile_im << x << " " << z << " " << psi[ir][it][ip].imag() << std::endl;
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
            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = theta.size() / 2; it < theta.size() / 2 + 1; ++it)
                {
                    for (std::size_t ip = 0; ip < phi.size(); ++ip)
                    {
                        number_t x = radial[ir] * sin(theta[it]) * cos(phi[ip]);
                        number_t y = radial[ir] * sin(theta[it]) * sin(phi[ip]);
                        myfile_phi_2 << x << " " << y << " " << psi[ir][it][ip].real() << std::endl;
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
            for (std::size_t ir = 0; ir < radial.size(); ++ir)
            {
                for (std::size_t it = theta.size() / 2; it < theta.size() / 2 + 1; ++it)
                {
                    for (std::size_t ip = 0; ip < phi.size(); ++ip)
                    {
                        number_t x = radial[ir] * sin(theta[it]) * cos(phi[ip]);
                        number_t y = radial[ir] * sin(theta[it]) * sin(phi[ip]);
                        myfile_phi_2_im << x << " " << y << " " << psi[ir][it][ip].imag() << std::endl;
                    }
                }
            }

            myfile_phi_2_im.close();
        }

        //void verify_solution(const Parameters& params, const QuantumState& qs)
        //{
        //    static constexpr complex_t CZERO{ 0.0, 0.0 };

        //    array3d<complex_t> ratio(radial.size(), phi.size(), theta.size());

        //    number_t alpha = params.h / (2.0 * params.mp);
        //    number_t E = params.ke / (qs.n * qs.n);

        //    for (std::size_t ir = 1; ir < radial.size() - 1; ++ir)
        //    {
        //        for (std::size_t it = 1; it < theta.size() - 1; ++it)
        //        {
        //            for (std::size_t ip = 1; ip < phi.size() - 1; ++ip)
        //            {
        //                if (psi[ir][it][ip] == CZERO)
        //                    continue;

        //                number_t dV1 = radial[ir] * radial[ir] * radial.delta;
        //                number_t dV2 = radial[ir] * radial[ir] * sin(theta[it]) * sin(theta[it]) * phi.delta;
        //                number_t dV3 = radial[ir] * radial[ir] * sin(theta[it]) * theta.delta;

        //                // number_t aS = r[i - 1] * r[i - 1] * alpha / (dV1 * r_p.delta);
        //                // number_t aN = r[i] * r[i] * alpha / (dV1 * r_p.delta);
        //                number_t r1 = radial[ir] - radial.delta / 2.0;
        //                number_t r2 = radial[ir] + radial.delta / 2.0;

        //                number_t aS = r1 * r1 * alpha / (dV1 * radial.delta);
        //                number_t aN = r2 * r2 * alpha / (dV1 * radial.delta);
        //                number_t aW = alpha / (dV2 * phi.delta);
        //                number_t aE = alpha / (dV2 * phi.delta);
        //                number_t aB = sin(theta[it] - theta.delta / 2.0) * alpha / (dV3 * theta.delta);
        //                number_t aT = sin(theta[it] + theta.delta / 2.0) * alpha / (dV3 * theta.delta);

        //                complex_t psi_s{ -aS * psi[ir - 1][ip][it] };
        //                complex_t psi_n{ -aN * psi[ir + 1][ip][it] };
        //                complex_t psi_w{ -aW * psi[ir][ip - 1][it] };
        //                complex_t psi_e{ -aE * psi[ir][ip + 1][it] };
        //                complex_t psi_b{ -aB * psi[ir][ip][it - 1] };
        //                complex_t psi_t{ -aT * psi[ir][ip][it + 1] };

        //                complex_t psi_p_k_rp{ -params.kp / radial[ir] * psi[ir][ip][it] };

        //                complex_t Fp{ psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + psi_p_k_rp };
        //                number_t ap = aS + aN + aW + aE + aB + aT - E;
        //                complex_t ap_psi_p{ -ap * psi[ir][it][ip] };

        //                ratio[ir][it][ip] = Fp / ap_psi_p;
        //            }
        //        }
        //    }

        //    /* Compute error */
        //    error = 0.0;
        //    for (std::size_t ir = 1; ir < radial.size() - 1; ++ir)
        //    {
        //        for (std::size_t it = 1; it < theta.size() - 1; ++it)
        //        {
        //            for (std::size_t ip = 1; ip < phi.size() - 1; ++ip)
        //            {
        //                if (fabs(ratio[ir][it][ip].real()) > 0.2)
        //                {
        //                    error += fabs(1.0 - ratio[ir][it][ip].real());
        //                }
        //            }
        //        }
        //    }
        //    error = 100.0 * error / ((radial.size() - 2) * (phi.size() - 2) * (theta.size() - 2));
        //}
    };
}

#endif /* !MANGH_HYDROGEN_SOLUTION_DATA_HPP */
