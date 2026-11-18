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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>

#include <CLI/CLI.hpp>

#include "parameters.hpp"
#include "quantum_state.hpp"
#include "mesh.hpp"
#include "math/constants.hpp"
#include "linspace.hpp"
#include "numbers.hpp"

// using namespace mangh::math;
using namespace mangh::hydrogenium;

/**
 * @brief Opens @p file_name for writing and passes the stream to @p export_func.
 *
 * @param file_name   - output file name,
 * @param export_func - any callable invocable as export_func(std::ostream&).
 */
template<typename ExportFunc>
void export_data(const char* file_name, ExportFunc&& export_func)
{
    std::ofstream myfile(file_name/* , std::ios::binary */);
    export_func(myfile);
    myfile.close();
    if (!myfile) std::cerr << "Failed to write '" << file_name << "'\n";
}

int main(int argc, char** argv)
{
    // get the program arguments
    CLI::App app{"Hydrogen Orbitals"};
    argv = app.ensure_utf8(argv);

    app.get_formatter()->column_width(35);
    app.get_formatter()->enable_option_type_names(false);
    app.get_formatter()->long_option_alignment_ratio( 0.2f /* 1/5*/);

    int principalQuantumNumber;
    app.add_option("-n,--principal", principalQuantumNumber, "Principal quantum number (1 <= n <= 10)")
        ->check(CLI::Range(1, 10))
        ->required();

    int angularQuantumNumber;
    app.add_option("-l,--angular", angularQuantumNumber, "Azimuthal number (0 <= l < n)")
        ->check(CLI::Range(0, 10))
        ->required();

    int magneticQuantumNumber;
    app.add_option("-m,--magnetic", magneticQuantumNumber, "Magnetic quantum number (-l <= m <= l)")
        ->check(CLI::Range(-10, 10))
        ->required();

    CLI11_PARSE(app, argc, argv);

    constexpr number_t PI = PI_V<number_t>;

    QuantumState qs
    {
        principalQuantumNumber,  // n - principal quantum number
        angularQuantumNumber,    // l - azimuthal number
        magneticQuantumNumber,   // m - magnetic quantum number
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
    export_data("grid_data.txt", [&](std::ostream& os) { mesh.export_grid_data(os /* , max_pd, max_real */); });

    /* Export probability density data */
    export_data("data.txt", [&](std::ostream& os) { mesh.export_probability_density(os); });

    /* Export psi real data */
    export_data("data_re.txt", [&](std::ostream& os) { mesh.export_psi_real(os); });

    /* Export psi im data */
    export_data("data_im.txt", [&](std::ostream& os) { mesh.export_psi_im(os); });

    /* Export psi real data at n_phi/2 */
    export_data("data_re_phi_2.txt", [&](std::ostream& os) { mesh.export_data_phi_2(os); });

    /* Export psi im data at n_phi/2 */
    export_data("data_im_phi_2.txt", [&](std::ostream& os) { mesh.export_data_phi_2_im(os); });

    ///* Verify solution numerically */
    //mesh.verify_solution(params, qs/* , radii, thetas, phis */);

    // char* prev_loc = std::setlocale(LC_ALL, nullptr);
    // std::setlocale(LC_ALL, "pl_PL.UTF-8");
    // if (const char* loc = std::setlocale(LC_ALL, "pl_PL.UTF-8"))
    //     std::printf("New LC_ALL locale: %s\n", loc);

    constexpr auto max_precision{std::numeric_limits<number_t>::digits10 + 1};     

    /* Print some results */
    for(std::size_t r = 0; r < mesh.radial.size(); ++r)
    {
        std::size_t theta = mesh.azimuthal.size() / 4;
        std::size_t phi = mesh.polar.size() / 2;
        std::cout
            << "real psi[" << r << "][" << theta << "][" << phi << "]: " << std::fixed << mesh.psi[r][theta][phi].real()
            << ", im psi[" << r << "][" << theta << "][" << phi << "]: " << std::fixed << mesh.psi[r][theta][phi].imag()
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
