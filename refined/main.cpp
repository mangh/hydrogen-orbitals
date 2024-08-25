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
    for(int i = 0; i < dd.n_r; ++i)
    {
        int theta = dd.n_theta/4;
        int phi = dd.n_phi/2;
        printf("real psi[%i][%i][%i]: %f, im psi[%i][%i][%i]: %f\n",
                i, theta, phi, sd.psi[i][theta][phi].real(),
                i, theta, phi, sd.psi[i][theta][phi].imag());
    }

    printf("error in percentage : %f\n", sd.error);
}

// #include <math.h>
// #include <stdio.h>
// // #include <cmath>
// #include "complex.hpp"
// #include "mem_ops.hpp"
// #include "support.hpp"
// #include "user_types.hpp"

// int main(/* int argc, char* argv[] */) {

//     d_data domain_data;
//     q_nums quantum_numbers;
//     p_params physical_params;
//     s_data solution_data;

//     /* Parameters */
//     domain_data.n_phi = 20;     //Number of nodes in phi direction
//     domain_data.n_theta = 40;   //Number of nodes in theta direction
//     domain_data.n_r = 200;      //Number of nodes in r direction
//     domain_data.R = 40.0;       //Length of domain

//     quantum_numbers.n = 4;      //Principal quantum number
//     quantum_numbers.l = 2;      //Azimuthal number
//     quantum_numbers.m = 1;      //Magnetic quantum number

//     physical_params.a0 = 1.0;   //Bohr radius
//     physical_params.h = 1.0;    //Constant
//     physical_params.mp = 1.0;   //Mass of particle
//     physical_params.kp = 1.0;   //Potential constant
//     physical_params.ke = 1.0;   //Constant in energy level equation

//     /* Allocate memory for solution data */
//     solution_data.psi = mat3D(domain_data.n_r, domain_data.n_theta, domain_data.n_phi);
//     solution_data.psi_square = mat3D(domain_data.n_r, domain_data.n_theta, domain_data.n_phi);
//     solution_data.r_p = new double[domain_data.n_r];
//     solution_data.theta_p = new double[domain_data.n_theta];
//     solution_data.phi_p = new double[domain_data.n_phi];

//     /* Compute wave function */
//     calc_psi(domain_data,
//              quantum_numbers,
//              physical_params,
//              &solution_data);

//     /* Print some results */
//     for(int i = 0; i < domain_data.n_r; ++i) {
//         int theta = domain_data.n_theta/4;
//         int phi = domain_data.n_phi/2;
//         printf("real psi[%i][%i][%i]: %f, im psi[%i][%i][%i]: %f\n",
//                 i, theta, phi, solution_data.psi[i][theta][phi].a,
//                 i, theta, phi, solution_data.psi[i][theta][phi].b);
//     }

//     printf("error in percentage : %f\n", solution_data.error);

//     /* Deallocate memory */
//     free_mat3D(solution_data.psi, domain_data.n_r, domain_data.n_theta);
//     free_mat3D(solution_data.psi_square, domain_data.n_r, domain_data.n_theta);
//     delete [] solution_data.r_p;
//     delete [] solution_data.theta_p;
//     delete [] solution_data.phi_p;

//     return 0;
// }
