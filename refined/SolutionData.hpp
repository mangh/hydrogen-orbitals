#ifndef SOLUTION_DATA_HPP
#define SOLUTION_DATA_HPP

// #include <complex>
#include <fstream>
// #include <vector>

#include "Data3D.hpp"
#include "DomainData.hpp"
#include "PhysicalParameters.hpp"
#include "QuantumNumbers.hpp"

namespace HydrogenOrbitals
{
    constexpr double M_PI = 3.141592653589793238462643383279502884L;

    // using ComplexData = std::complex<double>;

    /**
     * @brief Solution data.
     * 
     */
    struct SolutionData
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief brief description
        const DomainData& dd;

        /// @brief brief description
        std::vector<double> r;

        /// @brief brief description
        std::vector<double> r_p;

        /// @brief brief description
        std::vector<double> theta_p;

        /// @brief brief description
        std::vector<double> phi;

        /// @brief brief description
        std::vector<double> phi_p;

        // /// @brief brief std::complex<double> psi[d.n_r][d.n_theta][d.n_phi]
        // std::vector<std::vector<std::vector<ComplexData>>> psi;

        // /// @brief brief description
        // std::vector<std::vector<std::vector<ComplexData>>> psi_square;

        /// @brief brief description
        Data3D psi;

        /// @brief brief description
        Data3D psi_square;

        /// @brief brief description
        double E;

        /// @brief brief description
        double error;

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
            Initialize();

            // for(int i = 0; i < dd.n_r; ++i)
            // {
            //     psi[i] = std::vector<std::vector<ComplexData>>(dd.n_theta);
            //     psi_square[i] = std::vector<std::vector<ComplexData>>(dd.n_theta);
            //     for(int j = 0; j < dd.n_theta; ++j)
            //     {
            //         psi[i][j] = std::vector<ComplexData>(dd.n_phi);
            //         psi_square[i][j] = std::vector<ComplexData>(dd.n_phi);
            //     }
            // }
        }

        void Initialize(/* DomainData dd */)
        {
            double dr = dd.R / dd.n_r;
            double dtheta = 2.0 * M_PI / dd.n_theta;
            double dphi = M_PI / (dd.n_phi-1);

            // r & r_p vector 
            for(int i = 0; i < dd.n_r; ++i) {  r_p[i] = i * dr; r[i] = r_p[i] + 0.5 * dr; }

            // theta_p vector
            for(int i = 0; i < dd.n_theta; ++i) { theta_p[i] = i * dtheta + 0.5 * dtheta; }

            // phi vector
            for(int i = 0; i < dd.n_phi - 1; ++i) { phi[i] = i * dphi + 0.5 * dphi; }
            
            // phi_p vector
            for(int i = 0; i < dd.n_phi; ++i) { phi_p[i] = i * dphi; }
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        int factorial(int n)
        {
            return (n <= 0) ? 1 : n * factorial(n - 1);
        }

        long double P_l(long double x, int l)
        {
            long double P_l_x = 0.0L;
            for(int k = 0; k < l/2 + 1; ++k)
            {
                P_l_x += 1.0 / 
                    pow(2.0, l) * 
                    pow(-1.0, k) *
                    factorial(2 * l - 2 * k) * 
                    pow(x, l - 2 * k) /
                    (factorial(k) * factorial(l - k) * factorial(l - 2 * k));
            }
            return P_l_x;
        }

        long double dm_Pl_dxm(long double x, int m, int l)
        {
            long double dx = 0.0001L;
            return
                (m == 0) ? P_l(x, l) :
                (m == 1) ? (P_l(x + dx, l) - P_l(x, l)) / dx :
                (dm_Pl_dxm(x + dx, m - 1, l) - dm_Pl_dxm(x, m - 1, l)) / dx;
        }

        double legendre_polynomial(int m, int l, double x)
        {
            return(m >= 0) ?
                pow(-1.0, m) * pow(1.0 - x * x, 0.5 * m) * dm_Pl_dxm(x, m, l) :
                pow(-1.0, (-m)) * factorial(l - (-m)) / factorial(l + (-m)) * legendre_polynomial((-m), l, x);
        }

        long double laguerre_polynomial(long double x, int a, int k)
        {
            return
                (k == 0) ? 1.0 :
                (k == 1) ? 1.0 + a - x :
                ((2.0 * (k - 1) + 1.0 + a - x) * laguerre_polynomial(x, a, k - 1) - (k - 1 + a) * laguerre_polynomial(x, a, k - 2 )) / k;
        }

        void compute_wave_function(const PhysicalParameters& param, const QuantumNumbers& qnum)
        {
            for(int i = 0; i < 1; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        psi[i][j][k] = ComplexData{ 0.0, 0.0 };
                    }
                }
            }

            for(int i = 1; i < dd.n_r; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        double rho = 2.0 * r_p[i] / (qnum.n * param.a0);

                        double beta = 
                            exp(-rho / 2.0) *
                            pow(rho, qnum.l) *
                            laguerre_polynomial(rho, 2 * qnum.l + 1, qnum.n - qnum.l - 1) *
                            legendre_polynomial(qnum.m, qnum.l, cos(phi_p[k]));

                        ComplexData e_im_theta{ cos(qnum.m * theta_p[j]), sin(qnum.m * theta_p[j]) };

                        psi[i][j][k] = e_im_theta * beta;
                    }
                }
            }
        }

        void normalize_psi()
        {
            ComplexData integral{ 0.0, 0.0 };

            double dr = dd.R / dd.n_r;
            double dtheta = 2.0 * M_PI / dd.n_theta;
            double dphi = M_PI / (dd.n_phi-1);

            /* Top pole */
            for(int i = 1; i < dd.n_r; ++i)
            {
                double dV = 2.0 / 3.0 * M_PI * (1 - cos(dphi / 2)) * (r[i] * r[i] * r[i] - r[i-1] * r[i-1] * r[i-1]);
                ComplexData psi_sq = (psi[i][0][0] * conj(psi[i][0][0])) * dV;
                integral += psi_sq;
            }

            /* Central nodes */
            for(int i = 1; i < dd.n_r; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = 1; k < dd.n_phi - 1; ++k)
                    {
                        double dV = r_p[i] * r_p[i] * sin(phi_p[k]) * dphi * dtheta * dr;
                        ComplexData psi_sq = (psi[i][j][k] * conj(psi[i][j][k])) * dV;
                        integral += psi_sq;
                    }
                }
            }

            /* Bottom pole */
            for(int i = 1; i < dd.n_r; ++i)
            {
                double dV = 2.0 / 3.0 * M_PI * (1 - cos(dphi / 2.0)) * (r[i] * r[i] * r[i] - r[i-1] * r[i-1] * r[i-1]);
                ComplexData psi_sq = (psi[i][0][dd.n_phi - 1] * conj(psi[i][0][dd.n_phi - 1])) * dV;
                integral += psi_sq;
            }

            /* Normalize */
            for(int i = 1; i < dd.n_r; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        psi[i][j][k] /= sqrt(integral.real());
                    }
                }
            }
        }

        void compute_probability_density()
        {
            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        psi_square[i][j][k] = psi[i][j][k] * conj(psi[i][j][k]);
                    }
                }
            }
        }

        double compute_max_pd()
        {
            double max_pd = -1e+8;
            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        if(psi_square[i][j][k].real() > max_pd)
                        {
                            max_pd = psi_square[i][j][k].real();
                        }
                    }
                }
            }
            return max_pd;
        }

        double compute_max_psi_real()
        {
            double max_real = -1e+8;
            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        if(psi[i][j][k].real() > max_real)
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
            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = 0; j < 1; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        double x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        double z = r_p[i] * cos(phi_p[k]);
                        myfile_pd << x << " " << z << " " << psi_square[i][j][k].real() << std::endl;
                    }
                }
            }

            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = dd.n_theta/2; j < dd.n_theta/2 + 1; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        double x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        double z = r_p[i] * cos(phi_p[k]);
                        myfile_pd << x << " " << z << " " << psi_square[i][j][k].real() << std::endl;
                    }
                }
            }

            myfile_pd.close();
        }

        void export_psi_real(std::string file_name)
        {
            std::ofstream myfile_real;
            myfile_real.open(file_name);

            int theta_slice = 0; // Should be a value between 0 and n_theta / 2 - 2

            /* Export psi central nodes */
            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = theta_slice; j < theta_slice + 1; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        double x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        double z = r_p[i] * cos(phi_p[k]);
                        myfile_real << x << " " << z << " " << psi[i][j][k].real() << std::endl;
                    }
                }
            }

            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = dd.n_theta/2 + theta_slice; j < dd.n_theta/2 + theta_slice + 1; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
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

            int theta_slice = 0; // Should be a value between 0 and n_theta / 2 - 2

            /* Export psi central nodes */
            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = theta_slice; j < theta_slice + 1; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        double x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        double z = r_p[i] * cos(phi_p[k]);
                        myfile_im << x << " " << z << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = dd.n_theta / 2 + theta_slice; j < dd.n_theta / 2 + theta_slice + 1; ++j)
                {
                    for(int k = 0; k < dd.n_phi; ++k)
                    {
                        double x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        double z = r_p[i] * cos(phi_p[k]);
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
            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = dd.n_phi / 2; k < dd.n_phi / 2 + 1; ++k)
                    {
                        double x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        double y = r_p[i] * sin(phi_p[k]) * sin(theta_p[j]);
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
            for(int i = 0; i < dd.n_r; ++i)
            {
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    for(int k = dd.n_phi/2; k < dd.n_phi/2 + 1; ++k)
                    {
                        double x = r_p[i] * sin(phi_p[k]) * cos(theta_p[j]);
                        double y = r_p[i] * sin(phi_p[k]) * sin(theta_p[j]);
                        myfile_phi_2_im << x << " " << y << " " << psi[i][j][k].imag() << std::endl;
                    }
                }
            }

            myfile_phi_2_im.close();
        }

        void verify_solution(PhysicalParameters params, QuantumNumbers qnum)
        {
            Data3D ratio(dd);

            double dr = dd.R / dd.n_r;
            double dtheta = 2.0 * M_PI / dd.n_theta;
            double dphi = M_PI / (dd.n_phi - 1);
            double alpha = params.h / (2.0 * params.mp);
            double E = params.ke / (qnum.n * qnum.n);

            for(int i = 1; i < dd.n_r - 1; ++i)
            {
                for(int j = 1; j < dd.n_theta - 1; ++j)
                {
                    for(int k = 1; k < dd.n_phi - 1; ++k)
                    {
                        double dV1 = r_p[i] * r_p[i] * dr;
                        double dV2 = r_p[i] * r_p[i] * sin(phi_p[k]) * sin(phi_p[k]) * dtheta;
                        double dV3 = r_p[i] * r_p[i] * sin(phi_p[k]) * dphi;

                        double aS = r[i-1] * r[i-1] * alpha/ (dV1 * dr);
                        double aN = r[i] * r[i] * alpha / (dV1 * dr);
                        double aW = alpha / (dV2 * dtheta);
                        double aE = alpha / (dV2 * dtheta);
                        double aB = sin(phi[k-1]) *alpha / (dV3 * dphi);
                        double aT = sin(phi[k]) * alpha / (dV3 * dphi);

                        ComplexData psi_s{ -aS * psi[i-1][j][k] };
                        ComplexData psi_n{ -aN * psi[i+1][j][k] };
                        ComplexData psi_w{ -aW * psi[i][j-1][k] };
                        ComplexData psi_e{ -aE * psi[i][j+1][k] };
                        ComplexData psi_b{ -aB * psi[i][j][k-1] };
                        ComplexData psi_t{ -aT * psi[i][j][k+1] };
                        ComplexData psi_p_k_rp{ -params.kp / r_p[i] * psi[i][j][k] };

                        ComplexData Fp { psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + psi_p_k_rp };
                        double ap = aS + aN + aW + aE + aB + aT - E;
                        ComplexData ap_psi_p{ -ap * psi[i][j][k] };

                        ratio[i][j][k] = Fp / ap_psi_p;
                    }
                }
            }

            /* Compute error */
            error = 0.0;
            for(int i = 1; i < dd.n_r - 1; ++i)
            {
                for(int j = 1; j < dd.n_theta - 1; ++j)
                {
                    for(int k = 1; k < dd.n_phi - 1; ++k)
                    {
                        if(fabs(ratio[i][j][k].real()) > 0.2)
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
