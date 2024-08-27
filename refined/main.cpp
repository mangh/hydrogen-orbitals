/*
 * main.cpp
 *
 *  Created on: Jan 29, 2021
 *      Author: d-w-h
 */

#include "DomainData.hpp"
#include "QuantumNumbers.hpp"
#include "PhysicalParameters.hpp"
#include "SolutionData.hpp"
#include <ios>
#include <iostream>

using namespace HydrogenOrbitals;

int main()
{
    DomainData dd
    {
        200,    // number of nodes in r direction
        40,     // number of nodes in theta direction
        20,     // number of nodes in phi direction
        40.0    // length of domain
    };

    QuantumNumbers qn
    {
        4,  // n - Principal quantum number
        2,  // l - Azimuthal number
        1,  // m - Magnetic quantum number
    };

    PhysicalParameters params
    {
         1.0,   // a0_ - Bohr radius,
         1.0,   // h_  - Planck (???) constant,
         1.0,   // mp_ - Mass of particle,
         1.0,   // kp_ - Potential constant (???),
         1.0    // ke_ - Constant in energy level equation (???).
    };

    SolutionData sd{ dd };

    /* Compute psi */
    sd.compute_wave_function(params, qn);

    /* Normalize psi */
    sd.normalize_psi();

    /* Compute probability density */
    sd.compute_probability_density();

    /* Compute max probability density */
    double max_pd = sd.compute_max_pd();

    /* Compute max psi_real */
    double max_real = sd.compute_max_psi_real();

    /* Export domain data */
    sd.export_grid_data("grid_data.txt", max_pd, max_real);

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
    sd.verify_solution(params, qn);

    /* Print some results */
    for(std::size_t i = 0; i < dd.n_r; ++i)
    {
        std::size_t theta = dd.n_theta / 4;
        std::size_t phi = dd.n_phi / 2;
        std::cout
            << "real psi[" << i << "][" << theta << "][" << phi << "]: " << std::fixed << sd.psi[i][theta][phi].real()
            << ", im psi[" << i << "][" << theta << "][" << phi << "]: " << std::fixed << sd.psi[i][theta][phi].imag()
            << std::endl
            ;
    }

    std::cout << "error in percentage : " << std::fixed << sd.error << std::endl;
}
