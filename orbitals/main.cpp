/*
 * main.cpp
 *
 *  Created on: Jan 29, 2021
 *      Author: d-w-h
 */

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <iomanip>
#include <limits>

#include "parameters.hpp"
#include "quantum_state.hpp"
#include "mesh.hpp"
#include "math/constants.hpp"
#include "math/linspace.hpp"
#include "numbers.hpp"

// using namespace mangh::math;
using namespace mangh::hydrogenium;

int main()
{
    constexpr number_t PI = PI_V<number_t>;

    QuantumState qs
    {
        4,  // n - principal quantum number
        2,  // l - azimuthal number
        1,  // m - magnetic quantum number
    };

    Parameters params
    {
        1.0,   // h_  - Planck constant (???),
        1.0,   // mp_ - mass of particle,
        1.0,   // kp_ - potential constant (???),
        1.0    // ke_ - constant in energy level equation (???).
    };

    Mesh mesh
    {
        linspace<number_t>{ 0.0, 40.0, 200 },       //  radial coordinates (r)
        linspace<number_t>{ 0.0, PI, 20 },		    //  polar coordinates (θ)
        linspace<number_t>{ 0.0, 2.0 * PI, 40 }     //  azimuthal coordinates (ϕ)
    };

    // mesh.export_psi("psi.txt", params, qs, radii, thetas, phis);

    /* Compute psi */
    mesh.compute_wave_function(qs);

    /* Normalize psi */
    mesh.normalize_psi();

    // /* Compute probability density */
    // mesh.compute_probability_density();

    // /* Compute max probability density */
    // number_t max_pd = mesh.compute_max_pd();

    // /* Compute max psi_real */
    // number_t max_real = mesh.compute_max_psi_real();

    /* Export domain data */
    mesh.export_grid_data("grid_data.txt"/* , max_pd, max_real */);

    /* Export probability density data */
    mesh.export_probability_density("data.txt");

    /* Export psi real data */
    mesh.export_psi_real("data_real.txt");

    /* Export psi im data */
    mesh.export_psi_im("data_im.txt");

    /* Export psi real data at n_phi/2 */
    mesh.export_data_phi_2("data_real_phi_2.txt");

    /* Export psi im data at n_phi/2 */
    mesh.export_data_phi_2_im("data_im_phi_2.txt");

    ///* Verify solution numerically */
    //mesh.verify_solution(params, qs/* , radii, thetas, phis */);

    // char* prev_loc = std::setlocale(LC_ALL, nullptr);
    // std::setlocale(LC_ALL, "pl_PL.UTF-8");
    // if (const char* loc = std::setlocale(LC_ALL, "pl_PL.UTF-8"))
    //     std::printf("New LC_ALL locale: %s\n", loc);

    constexpr auto max_precision{std::numeric_limits<number_t>::digits10 + 1};     

    /* Print some results */
    for(std::size_t i = 0; i < mesh.radial.size(); ++i)
    {
        std::size_t theta = mesh.phi.size() / 4;
        std::size_t phi = mesh.theta.size() / 2;
        std::cout
            << "real psi[" << i << "][" << theta << "][" << phi << "]: " << std::fixed << mesh.psi[i][theta][phi].real()
            << ", im psi[" << i << "][" << theta << "][" << phi << "]: " << std::fixed << mesh.psi[i][theta][phi].imag()
            << std::endl
            ;
    }

    // // constexpr auto max_precision{std::numeric_limits<number_t>::digits10 + 1};     
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
    // //         << mesh.psi[i][theta][phi].real() << "\t"
    // //         << mesh.psi[i][theta][phi].imag()
    // //         << std::endl
    // //         ;
    // // }

    std::cout
        << std::scientific
        << std::setprecision(max_precision)
        << std::fixed
        //<< "error in percentage : " << mesh.error << std::endl
        << "normalization constant : " << std::sqrt(mesh.total_probability) << std::endl;

    // if (const char* loc = std::setlocale(LC_ALL, prev_loc/*.c_str()*/))
    //     std::printf("Restorred LC_ALL locale: %s\n", loc);

    return EXIT_SUCCESS;
}
