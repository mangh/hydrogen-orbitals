/*
 * main.cpp
 *
 *  Created on: Jan 29, 2021
 *      Author: d-w-h
 */

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

#include <CLI/CLI.hpp>

#include "quantum_state.hpp"
#include "mesh.hpp"
#include "numbers.hpp"

using namespace mangh::hydrogenium;

int main(int argc, char** argv)
try
{
    // get the program arguments
    CLI::App app{"Hydrogen Orbitals"};
    argv = app.ensure_utf8(argv);

    app.get_formatter()->column_width(35);
    app.get_formatter()->enable_option_type_names(false);
    app.get_formatter()->long_option_alignment_ratio( 0.2f /* 1/5*/);

    constexpr int max_n = QuantumState::max_n;

    int principalQuantumNumber;
    app.add_option("-n,--principal", principalQuantumNumber,
                   "Principal quantum number (1 <= n <= " + std::to_string(max_n) + ")")
        ->check(CLI::Range(1, max_n))
        ->required();

    int angularQuantumNumber;
    app.add_option("-l,--angular", angularQuantumNumber, "Azimuthal number (0 <= l < n)")
        ->check(CLI::Range(0, max_n - 1))
        ->required();

    int magneticQuantumNumber;
    app.add_option("-m,--magnetic", magneticQuantumNumber, "Magnetic quantum number (-l <= m <= l)")
        ->check(CLI::Range(-(max_n - 1), max_n - 1))
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

    const QuantumState qs
    {
        principalQuantumNumber,  // n - principal quantum number
        angularQuantumNumber,    // l - azimuthal number
        magneticQuantumNumber,   // m - magnetic quantum number
    };

    const Mesh mesh = Mesh::for_state(qs);

    /* Save the mesh (coordinates and Ψ) for visualization */
    std::filesystem::create_directories(outputDir);
    mesh.save(outputDir);
    std::cout << "Mesh saved to " << std::filesystem::absolute(outputDir).string() << std::endl;

    /* Print the total probability (≈ 1 for the normalized Ψ; the deficit is the quadrature error plus the tail beyond the grid) */
    constexpr auto max_precision{std::numeric_limits<real_t>::digits10 + 1};
    std::cout
        << std::setprecision(max_precision)
        << "Total probability over the grid: P = " << mesh.total_probability
        << std::setprecision(3) << std::scientific
        << "  (1 - P = " << 1.0 - mesh.total_probability << ")"
        << std::endl;

    return EXIT_SUCCESS;
}
catch (const std::exception& e)
{
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
}
