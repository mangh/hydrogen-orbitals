/*
 * main.cpp
 *
 *  Created on: Jan 29, 2021
 *      Author: d-w-h
 */

#include <ios>
#include <iostream>
#include <iomanip>

#include "PhysicalParameters.hpp"
#include "QuantumState.hpp"
#include "SolutionData.hpp"
#include "linspace.hpp"

using namespace HydrogenOrbitals;

int main()
{
    QuantumState qn
    {
        4,  // n - principal quantum number
        2,  // l - azimuthal number
        1,  // m - magnetic quantum number
    };

    PhysicalParameters params
    {
        1.0,   // a0_ - Bohr radius,
        1.0,   // h_  - Planck constant (???),
        1.0,   // mp_ - mass of particle,
        1.0,   // kp_ - potential constant (???),
        1.0    // ke_ - constant in energy level equation (???).
    };

    SolutionData sd
    {
		linspace<numeral_t>{ 0.0, 40.0, 200 },      //  radius
        linspace<numeral_t>{ 0.0, 2.0 * M_PI, 40 }, //  theta
		linspace<numeral_t>{ 0.0, M_PI, 20 }		//  phi
    };

    // sd.export_psi("psi.txt", params, qn, radii, thetas, phis);

    /* Compute psi */
    sd.compute_wave_function(params, qn);

    /* Normalize psi */
    sd.normalize_psi();

    // /* Compute probability density */
    // sd.compute_probability_density();

    // /* Compute max probability density */
    // numeral_t max_pd = sd.compute_max_pd();

    // /* Compute max psi_real */
    // numeral_t max_real = sd.compute_max_psi_real();

    /* Export domain data */
    sd.export_grid_data("grid_data.txt"/* , max_pd, max_real */);

    /* Export probability density data */
    sd.export_probability_density("data.txt");

    /* Export psi real data */
    sd.export_psi_real("data_real.txt");

    /* Export psi im data */
    sd.export_psi_im("data_im.txt");

    /* Export psi real data at n_phi/2 */
    sd.export_data_phi_2("data_real_phi_2.txt");

    /* Export psi im data at n_phi/2 */
    sd.export_data_phi_2_im("data_im_phi_2.txt");

    /* Verify solution numerically */
    sd.verify_solution(params, qn/* , radii, thetas, phis */);

    // char* prev_loc = std::setlocale(LC_ALL, nullptr);
    // std::setlocale(LC_ALL, "pl_PL.UTF-8");
    // if (const char* loc = std::setlocale(LC_ALL, "pl_PL.UTF-8"))
    //     std::printf("New LC_ALL locale: %s\n", loc);

    constexpr auto max_precision{std::numeric_limits<numeral_t>::digits10 + 1};     

    /* Print some results */
    for(std::size_t i = 0; i < sd.radial.size(); ++i)
    {
        std::size_t theta = sd.theta.size() / 4;
        std::size_t phi = sd.phi.size() / 2;
        std::cout
            << "real psi[" << i << "][" << theta << "][" << phi << "]: " << std::fixed << sd.psi[i][theta][phi].real()
            << ", im psi[" << i << "][" << theta << "][" << phi << "]: " << std::fixed << sd.psi[i][theta][phi].imag()
            << std::endl
            ;
    }

    // // constexpr auto max_precision{std::numeric_limits<numeral_t>::digits10 + 1};     
    // std::cout
    //     << "R = " << radii.back() << "\t"
    //     << "nr = " << dd.n_r << " (" << radii.size() << ")\t"
    //     << "n_theta = " << dd.n_theta << " (" << thetas.size() << ")\t"
    //     << "n_phi = " << dd.n_phi << " (" << phis.size() << ")"
    //     << std::endl
        
    //     << "r\ttheta\tphi\tpsi_re\tpsi_im"
    //     << std::endl
    //     ;

    // // https://stackoverflow.com/questions/15220861/how-can-i-set-the-decimal-separator-to-be-a-comma
    // std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));

    // // for(std::size_t i = 0; i < dd.n_r; ++i)
    // // {
    // //     std::size_t theta = dd.n_theta / 4;
    // //     std::size_t phi = dd.n_phi / 2;
    // //     std::cout
    // //         << i << "\t"
    // //         << theta << "\t"
    // //         << phi << "\t"
    // //         << std::scientific
    // //         << std::setprecision(max_precision)
    // //         << sd.psi[i][theta][phi].real() << "\t"
    // //         << sd.psi[i][theta][phi].imag()
    // //         << std::endl
    // //         ;
    // // }

    std::cout
        << std::scientific
        << std::setprecision(max_precision)
        << std::fixed
        << "error in percentage : " << sd.error << std::endl
        << "normalization constant : " << std::sqrt(sd.integral) << std::endl;

    // if (const char* loc = std::setlocale(LC_ALL, prev_loc/*.c_str()*/))
    //     std::printf("Restorred LC_ALL locale: %s\n", loc);
}
