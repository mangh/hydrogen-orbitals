/*
 * main.cpp
 *
 *  Created on: Jan 29, 2021
 *      Author: d-w-h
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ios>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

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
    app.add_option("-n,--principal", principalQuantumNumber, "Principal quantum number (1 <= n <= 20)")
        ->check(CLI::Range(1, 20))
        ->required();

    int angularQuantumNumber;
    app.add_option("-l,--angular", angularQuantumNumber, "Azimuthal number (0 <= l < n)")
        ->check(CLI::Range(0, 19))
        ->required();

    int magneticQuantumNumber;
    app.add_option("-m,--magnetic", magneticQuantumNumber, "Magnetic quantum number (-l <= m <= l)")
        ->check(CLI::Range(-19, 19))
        ->required();

    // cross-option checks: run once after all options are parsed and stored
    app.final_callback([&]() {
        if (angularQuantumNumber >= principalQuantumNumber)
            throw CLI::ValidationError("--angular",
                "l = " + std::to_string(angularQuantumNumber) +
                " must satisfy 0 <= l < n (n = " + std::to_string(principalQuantumNumber) + ")");
        if (std::abs(magneticQuantumNumber) > angularQuantumNumber)
            throw CLI::ValidationError("--magnetic",
                "m = " + std::to_string(magneticQuantumNumber) +
                " must satisfy -l <= m <= l (l = " + std::to_string(angularQuantumNumber) + ")");
    });

    CLI11_PARSE(app, argc, argv);

    constexpr real_t PI = PI_V<real_t>;

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

    // Radial range scales with n (see QuantumState::RadialExtent), with a fixed step
    constexpr real_t radialStep = 0.2;
    const real_t radialMax = qs.RadialExtent();
    const auto radialPoints = static_cast<std::size_t>(std::lround(radialMax / radialStep)) + 1;

    // Polar grid θ ∈ [0, π]:
    // - an odd number of points, so that the slice index polarPoints/2 lies exactly at θ = π/2 (XY-plane exports),
    // - N = polarPoints - 1 >= 2l intervals: the θ integral (Clenshaw–Curtis, see Mesh::compute_total_probability) is then exact for |Yl,m|².
    const std::size_t polarPoints = std::max<std::size_t>(21, 2 * qs.l + 1);

    // Azimuthal grid ϕ ∈ [0, 2π] (the last node duplicates the first):
    // - an even number of intervals, so that the index azimuthalPoints/2 lies exactly at ϕ = π (opposite half-plane in the XZ-plane exports),
    // - >= 4|m| intervals: at least 4 samples per period of cos(mϕ) in the XY-plane plots
    //   (and > 2|m|, which also keeps the ϕ integral exact for real orbitals, where |Ψ|² contains cos²(mϕ)).
    const std::size_t azimuthalPoints = std::max<std::size_t>(40, static_cast<std::size_t>(4 * std::abs(qs.m))) + 1;

    Mesh mesh
    {
        linspace<real_t>{ 0.0, radialMax, radialPoints },       //  radial coordinates (r)
        linspace<real_t>{ 0.0, PI, polarPoints },               //  polar coordinates (θ)
        linspace<real_t>{ 0.0, 2.0 * PI, azimuthalPoints },  //  azimuthal coordinates (ϕ)
        qs
    };

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

    // char* prev_loc = std::setlocale(LC_ALL, nullptr);
    // std::setlocale(LC_ALL, "pl_PL.UTF-8");
    // if (const char* loc = std::setlocale(LC_ALL, "pl_PL.UTF-8"))
    //     std::printf("New LC_ALL locale: %s\n", loc);

    constexpr auto max_precision{std::numeric_limits<real_t>::digits10 + 1};     

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

    // // constexpr auto max_precision{std::numeric_limits<real_t>::digits10 + 1};     
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
        << "normalizing constant : " << std::sqrt(mesh.normalizing_factor) << std::endl;

    // if (const char* loc = std::setlocale(LC_ALL, prev_loc/*.c_str()*/))
    //     std::printf("Restorred LC_ALL locale: %s\n", loc);

    return EXIT_SUCCESS;
}
