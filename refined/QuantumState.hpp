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

        static bool permissible(unsigned n, unsigned l, int m, double s)
        {
            return (0 < n) && (l < n) && (static_cast<unsigned>(abs(m)) <= l) && ((s == HALFSPIN) || (s == -HALFSPIN));
        }

        static void check(unsigned n, unsigned l, int m, double s)
        {
            if (!permissible(n, l, m, s))
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
         * @brief Wave function Ψ (non-normalized).
         *
         * @param rho - radius ρ = 2 * r / (n * a0),
         * @param theta - azimuthal angle θ (the angle in the xy-plane from the x-axis with 0 ≤ θ < 2π),
         * @param phi - polar angle ϕ (the angle from the positive z-axis with 0 ≤ ϕ ≤ π),
         *
         * @remarks Ψ[n,l,m](ρ,θ,ϕ) = C * e^(-ρ/2) * ρ^l * L[n-l-1,2l+1](ρ) * e^imθ * P[l,m](cos⁡ ϕ)
         *
         */
        complex_t PSI(numeral_t rho, numeral_t theta, numeral_t phi) const
        {
            // return exp(-rho / 2.0)
            //      * pow(rho, l)
            //      * pow(-1, m) * std::assoc_laguerre(n - l - 1, 2 * l + 1, rho)
            //      * std::assoc_legendre(l, m, cos(phi))
            //      * complex_t{ cos(m * theta), sin(m * theta) };
            return exp(-rho / 2.0)
                 * pow(rho, l)
                 * assoc_laguerre(n - l - 1, 2 * l + 1, rho)
                 * assoc_legendre(l, m, cos(phi))
                 * complex_t{ cos(m * theta), sin(m * theta) };
        }
    };
}

#endif /* !QUANTUM_NUMBERS_HPP */