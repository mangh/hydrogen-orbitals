/*
 * main.cpp
 *
 *  Created on: Jan 29, 2021
 *      Author: d-w-h
 */

#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <iomanip>
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

    std::string basisName = "complex";
    app.add_option("-b,--basis", basisName,
                   "Angular functions: complex (e^(im*phi), eigenstates of Lz) or real (px, py, dxy, ...)")
        ->check(CLI::IsMember({ "complex", "real" }))
        ->capture_default_str();

    Resolution resolution;
    app.add_option("--radial-step", resolution.radial_step, "Radial step in Bohr radii")
        ->default_str("0.2");
    app.add_option("--polar-points", resolution.polar_points,
                   "Polar nodes theta in [0, pi]: odd, >= 2l+1 (default: max(21, 2l+1))");
    app.add_option("--azimuthal-points", resolution.azimuthal_points,
                   "Azimuthal nodes phi in [0, 2pi] incl. 2pi: odd, >= 5, > 2|m|+1 (default: max(40, 4|m|)+1)");

    std::string outputDir;
    app.add_option("-o,--output", outputDir, "Output directory (default: results/n<n>_l<l>_m<m>[_real])");

    CLI11_PARSE(app, argc, argv);

    const Basis basis = basisName == "real" ? Basis::Real : Basis::Complex;

    if (outputDir.empty())
        outputDir = "results/n" + std::to_string(principalQuantumNumber)
                  + "_l" + std::to_string(angularQuantumNumber)
                  + "_m" + std::to_string(magneticQuantumNumber)
                  + (basis == Basis::Real ? "_real" : "");

    const QuantumState qs
    {
        principalQuantumNumber,  // n - principal quantum number
        angularQuantumNumber,    // l - azimuthal number
        magneticQuantumNumber,   // m - magnetic quantum number
        basis                    // complex or real spherical harmonics
    };

    const Mesh mesh = Mesh::for_state(qs, resolution);

    /* Save the mesh (coordinates and Ψ) for visualization */
    std::filesystem::create_directories(outputDir);
    mesh.save(outputDir);

    // Console output in ASCII: the Windows console may not be able to print Greek letters.
    std::cout
        << "Orbital " << qs.Name() << ": n = " << qs.n << ", l = " << qs.l << ", m = " << qs.m
        << ", " << basisName << " basis\n"
        << "Mesh: " << mesh.psi.pages() << " x " << mesh.psi.rows() << " x " << mesh.psi.columns()
        << " nodes (r, theta, phi), r <= " << mesh.radial.back() << ", radial step " << mesh.radial.delta
        << ", psi.npy " << std::fixed << std::setprecision(1)
        << static_cast<double>(mesh.psi.size() * sizeof(complex_t)) / 1e6 << " MB\n"
        << "Saved to " << std::filesystem::absolute(outputDir).string() << "\n\n";

    /* Checks: the total probability (≈ 1 for the normalized Ψ; the deficit is the quadrature error plus the tail
       beyond the grid), and expectation values and energy on the mesh against the exact formulas */
    auto row = [](const char* label, real_t grid, real_t exact)
    {
        std::cout << "  " << std::left << std::setw(26) << label << std::right
                  << std::defaultfloat << std::setprecision(12) << std::setw(20) << grid
                  << std::setw(20) << exact
                  << std::scientific << std::setprecision(2) << std::setw(14)
                  << std::abs(grid - exact) / std::abs(exact) << "\n";
    };

    constexpr real_t hartree_eV = 27.211386245988;  // CODATA 2018
    const real_t energy = -mesh.mean_inv_r / 2.0;    // virial theorem: <T> = -E, <V> = 2E = -<1/r>

    std::cout << "Checks" << std::setw(43) << "on the mesh" << std::setw(20) << "exact" << std::setw(14) << "rel. error" << "\n";
    row("total probability", mesh.total_probability, 1.0);
    row("<r>", mesh.mean_r, qs.MeanR());
    row("<r^2>", mesh.mean_r2, qs.MeanR2());
    row("<1/r>", mesh.mean_inv_r, qs.MeanInvR());
    row("energy -<1/r>/2 [Eh]", energy, qs.Energy());
    row("energy -<1/r>/2 [eV]", energy * hartree_eV, qs.Energy() * hartree_eV);
    std::cout << std::flush;

    return EXIT_SUCCESS;
}
catch (const std::exception& e)
{
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
}
