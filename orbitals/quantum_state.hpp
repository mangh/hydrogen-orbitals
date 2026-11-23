#ifndef MANGH_HYDROGENIUM_STATE_HPP
#define MANGH_HYDROGENIUM_STATE_HPP

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "math/assoc_laguerre.hpp"
#include "math/factorial.hpp"
#include "math/spherical_harmonic.hpp"

#include "numbers.hpp"

namespace mangh { namespace hydrogenium
{
    using namespace mangh::math;
    using std::invalid_argument;

    /**
     * @brief Basis of the angular functions.
     *
     * - `Complex` - the complex spherical harmonics `Yl,m(θ,ϕ) ∝ e^(imϕ)`: eigenfunctions of Lz (m = the Lz eigenvalue),
     *
     * - `Real` - the real spherical harmonics (the orbitals of chemistry: px, py, pz, dxy, ...):
     *   `∝ cos(mϕ)` for m > 0, `∝ sin(|m|ϕ)` for m < 0, and `Yl,0` for m = 0.
     *   They are real combinations of `Yl,m` and `Yl,-m`, so they are not eigenfunctions of Lz (except for m = 0).
     */
    enum class Basis { Complex, Real };

    /**
     * @brief Quantum state.
     * 
     */
    struct QuantumState
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Physical Constants
        ///

        // h    6.62607015e-34         [J/Hz]  Planck const
        // ℏ    1.054571817e−34        [J*s]   Planck const reduced ℏ = h/2π
        // c    299792458              [m/s]   Speed of light
        // ε0   8.8541878188e-12       [F/m]   Vacuum permittivity constant

        // me   9.1093837139e-31       [kg]    Mass of electron
        // mp   1.67262192595e-27      [kg]    Mass of proton
        // e    1.602176634e-19        [C]     Elementary charge

        // α    0.0072973525643(11)    []      Fine structure constant = e^2 / (4·π·ε0·ℏ·c)
        // a0   5.291772105443463e-11  [m]     Bohr radius = 4*π*ε0*ℏ^2 / (me*e^2)
        // μ    9.10442528891677e-31   [kg]    Reduced mass of the electron–proton system = me * mp / (me + mp)
        // ra0  5.294654094605935e-11  [m]     Bohr radius reduced for the reduced mass of the electron-proton system = (me / μ) * a0

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Limits
        ///

        /// @brief The largest principal quantum number supported
        /// (RadialExtent() was verified up to this n, see its remarks).
        static constexpr int max_n = 20;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        QuantumState() = delete;

        /**
         * @brief QuantumState constructor.
         * 
         * @param n - principal quantum number (n > 0),
         * @param l - azimuthal number (0 <= l < n),
         * @param m - magnetic quantum number (-l <= m <= l),
         * @param basis - basis of the angular functions (complex or real spherical harmonics).
         *
         * @exception std::invalid_argument - raised when arguments (n, l, m) do not meet the prerequisites.
         * 
         */
        QuantumState(int n, int l, int m, Basis basis = Basis::Complex)
            : n(n > 0 ? static_cast<unsigned>(n) : throw std::invalid_argument(badQuantumNumbers))
            , l(0 <= l && l < n ? static_cast<unsigned>(l) : throw std::invalid_argument(badQuantumNumbers))
            , m(-l <= m && m <= l ? m : throw std::invalid_argument(badQuantumNumbers))
            , basis(basis)
            // auxiliary variables
            , ra0(1.0)
            , nra0(ra0 * n)
            , radial_factor(
                std::sqrt(
                    std::pow(2.0 / nra0, 3) * factorial<real_t>(static_cast<unsigned>(n - l - 1)) / 2.0 / n / factorial<real_t>(static_cast<unsigned>(n + l))
                )
            )
        {
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        /**
         * @brief The normalized radial part `Rn,l(r)` of the wave function `Ψ`.
         *
         * - `Rn,l(r) = N · e^(-ρ) · (2ρ)^l · L_n-l-1_2l+1(2ρ)`,
         *
         * - `ρ = r / nra0 = r / (n·a0)` - scaled radial distance from the nucleus,
         *
         * - `N = √[ (2/(n·a0))³ · (n-l-1)! / (2n·(n+l)!) ]` - normalizing factor.
         *
         * @param r - radial distance from the nucleus.
         *
         * @returns The value of the function `Rn,l` at `(r)`.
         * 
         */
        real_t Radial(real_t r) const
        {
            auto rho = r / nra0;
            return
                radial_factor * 
                std::exp(-rho) * 
                std::pow(2 * rho, l) *
                assoc_laguerre(2 * l + 1, n - l - 1, 2 * rho);
        }

        /**
         * @brief Radial extent of the orbital: the radius enclosing all but < 5e-6 of the probability for any l (0 <= l < n, n <= max_n).
         *
         * - `R_max = ra0 · n · (2n + 12)`
         *
         * @remarks The 5e-6 bound was verified numerically for n = 1..max_n
         * (the tail outside `R_max` slowly grows with n, so re-check it before raising max_n;
         * the test `radial_extent` in tests/test_orbitals.cpp does that).
         *
         * - `2n²` is the outermost classical turning point (for l = 0),
         * - `12n` is a margin proportional to the length of the probability density decay (which scales with n).
         *
         * @returns The radius beyond which the wave function can be neglected.
         *
         */
        real_t RadialExtent() const
        {
            return nra0 * (2.0 * n + 12.0);
        }

        /**
         * @brief Angular part of the wave function (aka spherical harmonics `Yl,m(θ,ϕ)`) in the state's basis.
         *
         * - complex basis: `Yl,m(θ,ϕ) = Pl,m(cosθ) · e^(imϕ)`,
         *
         * - real basis: AngularRealBasis(θ,ϕ) (returned as a complex number with a zero imaginary part).
         *
         * @param theta - polar angle 0 ≤ θ ≤ π (from the positive Z axis),
         * @param phi - azimuthal angle 0 ≤ ϕ < 2π (from the X-axis, in the XY-plane).
         *
         * @returns The value of the angular function `Yl,m` at `(θ,ϕ)`.
         * 
         * @remarks The calculation uses the relationship (6.7.3):
         * 
         * - `Yl,-m(θ,ϕ) = (-1)^m · conj(Yl,m(θ,ϕ))`
         *
         * from Numerical Recipes in C++, Chapter 6.7, Spherical Harmonics.
         * 
         */
        complex_t Angular(real_t theta, real_t phi) const
        {
            return basis == Basis::Real
                ? complex_t(spherical_real_basis(l, m, theta, phi))
                : spherical_harmonic(l, m, theta, phi);
        }

        /**
         * @brief Real part of the angular function (spherical harmonics in real basis).
         *
         * @param theta - polar angle 0 ≤ θ ≤ π (from the positive Z axis),
         * @param phi - azimuthal angle 0 ≤ ϕ < 2π (from the X-axis, in the XY-plane).
         *
         * @returns The value of the real angular function `Yl,m` at `(θ,ϕ)`.
         *
         */
        real_t AngularRealBasis(real_t theta, real_t phi) const
        {
            return spherical_real_basis(l, m, theta, phi);
        }
 
        /**
         * @brief Normalized wave function.
         *
         * - `Ψn,l,m(r,θ,ϕ) = Rn,l(r) * Yl,m(θ,ϕ)`.
         *
         * @param r - radial distance from the nucleus,
         * @param theta - polar angle 0 ≤ θ ≤ π (from the positive Z axis),
         * @param phi - azimuthal angle 0 ≤ ϕ < 2π (from the X-axis, in the XY-plane).
         *
         * @returns The value of the wave function (probability amplitude).
         *
         */
        complex_t WaveFunction(real_t r, real_t theta, real_t phi) const
        {
            return Radial(r) * Angular(theta, phi);
        }

        /**
         * @brief Name of the orbital (ASCII), e.g. `1s`, `2p_x`, `3d_z2`, `4f_xyz` (real basis), `4f(m=-1)` (complex basis).
         *
         * @remarks The real orbitals of l <= 3 are named after the Cartesian polynomials they are proportional to
         * (`z2` stands for 3z² - r², `x2-y2` for x² - y², and so on); for l > 3 the name gives m (`5g_m-2`).
         * The letters for l are s, p, d, f, g, h, i, k, ... (alphabetical from f on, without j, and without the p and s already used).
         */
        std::string Name() const
        {
            static constexpr const char* letters = "spdfghiklmnoqrtuvwxy";
            static constexpr const char* cartesian[4][7] = {
                { "", "", "", "", "", "", "" },
                { "y", "z", "x", "", "", "", "" },
                { "xy", "yz", "z2", "xz", "x2-y2", "", "" },
                { "y(3x2-y2)", "xyz", "yz2", "z3", "xz2", "z(x2-y2)", "x(x2-3y2)" },
            };

            std::string name = std::to_string(n) + letters[l];
            if (l == 0)
                return name;
            if (basis == Basis::Complex)
                return name + "(m=" + std::to_string(m) + ")";
            if (l <= 3)
                return name + "_" + cartesian[l][m + static_cast<int>(l)];
            return name + "_m" + std::to_string(m);
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Exact values (to check the numerical results against)
        ///
        ///     The radial expectation values are the same in both bases (they depend on n and l only).
        ///

        /// @brief Energy `En = -1/(2n²)` in hartrees (`1 Eh = 27.211386245988 eV`), with a0 = ra0 = 1.
        real_t Energy() const { return -1.0 / (2.0 * n * n * ra0); }

        /// @brief Mean distance from the nucleus `<r> = a0 · (3n² - l(l+1)) / 2`.
        real_t MeanR() const { return ra0 * (3.0 * n * n - l * (l + 1.0)) / 2.0; }

        /// @brief Mean squared distance `<r²> = a0² · n² (5n² + 1 - 3l(l+1)) / 2`.
        real_t MeanR2() const { return ra0 * ra0 * n * n * (5.0 * n * n + 1.0 - 3.0 * l * (l + 1.0)) / 2.0; }

        /// @brief Mean inverse distance `<1/r> = 1 / (a0 · n²)`; by the virial theorem `E = -<1/r> / 2`.
        real_t MeanInvR() const { return 1.0 / (ra0 * n * n); }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Principal quantum number n (n > 0).
        /// Represents the electron's energy level and relative size of the orbital.
        const unsigned n;

        /// @brief Azimuthal quantum number l (0 <= l < n).
        /// Relates to the shape of the atomic orbital.
        const unsigned l;

        /// @brief Magnetic quantum number m (-l <= m <= l).
        /// Specifies the orientation of the orbital in space.
        const int m;

        /// @brief Basis of the angular functions (complex or real spherical harmonics).
        const Basis basis;
 
    private:

        /// @brief Reduced Bohr radius
        /// for the reduced mass of the electron-proton system.
        const real_t ra0;

        /// @brief n * reduced Bohr radius (= n * ra0).
        const real_t nra0;

        /// @brief Normalizing factor
        /// for the radial part Rn,l of the wave function.
        const real_t radial_factor;

        static constexpr const char * badQuantumNumbers =
            "Quantum numbers (n, l, m) do not meet the conditions: n > 0,  0 <= l < n, -l <= m <= l.";
    };
}}

#endif /* !MANGH_HYDROGENIUM_STATE_HPP */
