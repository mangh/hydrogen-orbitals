#ifndef QUANTUM_NUMBERS_HPP
#define QUANTUM_NUMBERS_HPP

#include <stdexcept>
#include <string>
#include <vector>

namespace HydrogenOrbitals
{
    using std::invalid_argument;
    using std::to_string;

    /**
     * @brief Quantum numbers.
     * 
     */
    struct QuantumNumbers
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
        /// @pre: 1 <= n.
        const int n;

        /// @brief Azimuthal quantum number (relates to the shape of the atomic orbital).
        /// @pre: 0 <= l <= n-1.
        const int l;

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
        QuantumNumbers() = delete;

        /**
         * @brief QuantumNumbers constructor.
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
        QuantumNumbers(int n, int l, int m, double s = HALFSPIN)
            : n(n)
            , l(l)
            , m(m)
            , s(s)
        {
            check(n, l, m, s);
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Validation methods
        ///

        static bool permissible(int n, int l, int m, double s)
        {
            return (1 <= n) &&
                   (0 <= l) && (l <= n - 1) &&
                   (-l <= m) && (m <= l) &&
                   ((s == HALFSPIN) || (s == -HALFSPIN));
        }

        static void check(int n, int l, int m, double s)
        {
            if (!permissible(n, l, m, s))
            {
                const char* format = "Invalid quantum numbers :: n = %d :: l = %d :: m = %d :: s = %f";
                std::size_t len = std::snprintf(nullptr, 0, format, n, l, m, s);
                std::vector<char> buf(len + 1);
                std::snprintf(&buf[0], buf.size(), format, n, l, m, s);
                throw invalid_argument(std::string(buf.cbegin(), buf.cend() - 1));
            }
        }
    };
}

#endif /* !QUANTUM_NUMBERS_HPP */