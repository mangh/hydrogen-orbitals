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
#include <filesystem>
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

    std::string outputDir;
    app.add_option("-o,--output", outputDir, "Output directory (default: results/n<n>_l<l>_m<m>)");

    CLI11_PARSE(app, argc, argv);

    if (outputDir.empty())
        outputDir = "results/n" + std::to_string(principalQuantumNumber)
                  + "_l" + std::to_string(angularQuantumNumber)
                  + "_m" + std::to_string(magneticQuantumNumber);

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
    // - an odd number of points, so that a node lies exactly at θ = π/2 (XY-plane plots),
    // - N = polarPoints - 1 >= 2l intervals: the θ integral (Clenshaw–Curtis, see Mesh::compute_total_probability) is then exact for |Yl,m|².
    const std::size_t polarPoints = std::max<std::size_t>(21, 2 * qs.l + 1);

    // Azimuthal grid ϕ ∈ [0, 2π] (the last node duplicates the first):
    // - an even number of intervals, so that for every node ϕ there is also a node at ϕ + π (the opposite half-plane in the XZ-plane plots),
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

    /* Save the mesh (coordinates and Ψ) for visualization */
    std::filesystem::create_directories(outputDir);
    mesh.save(outputDir);
    std::cout << "Mesh saved to " << std::filesystem::absolute(outputDir).string() << std::endl;

    /* Print the normalization factor */
    constexpr auto max_precision{std::numeric_limits<real_t>::digits10 + 1};     
    std::cout
        << std::scientific
        << std::setprecision(max_precision)
        << std::fixed
        << "Normalization factor, i.e., the total probability for the entire volume defined by the grid: "
        << std::sqrt(mesh.normalizing_factor)
        << std::endl;

    return EXIT_SUCCESS;
}
