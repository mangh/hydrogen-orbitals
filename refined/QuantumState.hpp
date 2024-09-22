#ifndef QUANTUM_NUMBERS_HPP
#define QUANTUM_NUMBERS_HPP

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "Numeral.hpp"
#include "SpecialPolynomial.hpp"

namespace HydrogenOrbitals
{
    using std::invalid_argument;
    using std::to_string;
    using namespace SpecialPolynomial;

    /**
     * @brief Quantum numbers.
     * 
     */
    struct QuantumState
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constants
        ///

        static constexpr double HALFSPIN = 1.0 / 2.0;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Principal quantum number (represents the electron's energy level and relative size of the orbital).
        /// @pre: 0 < n.
        const unsigned n;

        /// @brief Azimuthal quantum number (relates to the shape of the atomic orbital).
        /// @pre: 0 <= l < n.
        const unsigned l;

        /// @brief Magnetic quantum number (specifies the orientation of the orbital in space).
        /// @pre: -l <= m <= l.
        const int m;

        /// @brief Spin quantum number (describes the electron's intrinsic spin).
        /// @remarks In a hydrogen-like atom, the spin of an electron does not affect the spatial distribution of the electron cloud.
        /// @pre: s == +1/2 || s == -1/2.
        const double s;

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
         * @param n - principal quantum number,
         * @param l - azimuthal number,
         * @param m - magnetic quantum number.
         * @param s - spin quantum number.
         * @exception invalid_argument exception description
         *
         * @remarks Data used in Derek's (original) version: n = 4, l = 2, m = 1.
         *
         */
        QuantumState(unsigned n, unsigned l, int m, double s = HALFSPIN)
            : n(n)
            , l(l)
            , m(m)
            , s(s)
        {
            check(n, l, m, s);
        }

        static bool is_permissible(unsigned n, unsigned l, int m, double s)
        {
            return (0 < n) && (l < n) && (static_cast<unsigned>(abs(m)) <= l) && ((s == HALFSPIN) || (s == -HALFSPIN));
        }

        static void check(unsigned n, unsigned l, int m, double s)
        {
            if (!is_permissible(n, l, m, s))
            {
                const char* format = "Invalid quantum numbers :: n = %u :: l = %u :: m = %d :: s = %f";
                std::size_t len = std::snprintf(nullptr, 0, format, n, l, m, s);
                std::vector<char> buf(len + 1);
                std::snprintf(&buf[0], buf.size(), format, n, l, m, s);
                throw invalid_argument(std::string(buf.cbegin(), buf.cend() - 1));
            }
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///

        /**
         * @brief Normalizing factor for the radial part (Rnl) of the wave function in the (n,l) state.
         * 
         * @param n - principal quantum number,
         * @param l - azimutha quantum number,
         * @param a0 - Bohr radius
         * @returns numeral_t 
         */
        static numeral_t Rnl_normalizing_factor(unsigned n, unsigned l, numeral_t a0)
        {
            return std::sqrt(
                std::pow(2.0 / n / a0, 3) *
                fact<numeral_t>(n - l - 1) /
                2.0 / n / fact<numeral_t>(n + l)
            );
        }

        /**
         * @brief Radial part Rn,l(r) of the wave function (not ortho-normalized).
         * @details Rn,l(ρ) = e^(-ρ) * (2ρ)^l * L[n-l-1,2l+1](2ρ).
         * 
         * @param rho - normalized radial distance from the nucleus: ρ = r / (n * a0).
         * @param rho - normalized radial distance from the nucleus: ρ = r / a0.
         * @returns The value of the radial function Rn,l(ρ).
         */
        numeral_t Rnl(numeral_t rho) const
        {
            return
                exp(-rho / n) *
                pow(2 * rho / n, l) * 
                assoc_laguerre(n - l - 1, 2 * l + 1, 2 * rho / n);
        }

        /**
         * @brief Angular part (spherical harmonics) Yl,m(θ,φ) of the wave function.
         * @details Yl,m(θ,φ) = Pl,m(cos⁡ ϕ) * e^imθ.
         * 
         * @param theta - azimuthal angle θ (the angle in the xy-plane from the x-axis with 0 ≤ θ < 2π),
         * @param phi - polar angle ϕ (the angle from the positive z-axis with 0 ≤ ϕ ≤ π),
         * @returns The value of the spherical harmonic Yl,m(θ,φ).
         *
         */
        complex_t Ylm(numeral_t theta, numeral_t phi) const
        {
            unsigned mm = (m < 0) ? -m : m;

            complex_t ylm = ylm_assoc_legendre(l, mm, cos(phi)) *
                complex_t{ cos(mm * theta), sin(mm * theta) };

            return (m >= 0) ? ylm : (mm & 1) ? -conj(ylm) : conj(ylm);
        }

        /**
         * @brief Wave function Ψn,l,m (non-normalized).
         * @details Ψn,l,m(ρ,θ,ϕ) = Rn,l(ρ) * Yl,m(θ, ϕ).
         *
         * @param rho - normalized radial distance from the nucleus: ρ = r / (n * a0).
         * @param theta - azimuthal angle θ (the angle in the xy-plane from the x-axis with 0 ≤ θ < 2π),
         * @param phi - polar angle ϕ (the angle from the positive z-axis with 0 ≤ ϕ ≤ π),
         *
         * @returns The value of the wave function Ψn,l,m(ρ,θ,ϕ).
         *
         */
        complex_t WaveFunction(numeral_t rho, numeral_t theta, numeral_t phi) const
        {
            return Rnl_normalizing_factor(n, l, 1.0) * Rnl(rho) * Ylm(theta, phi);
        }
    };
}

#endif /* !QUANTUM_NUMBERS_HPP */