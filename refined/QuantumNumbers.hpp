#ifndef QUANTUM_NUMBERS_HPP
#define QUANTUM_NUMBERS_HPP

#include <stdexcept>
#include <string>

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
         *
         * @remarks Data used in Derek's (original) version: n = 4, l = 2, m = 1.
         *
         */
        QuantumNumbers(int n, int l, int m, double s = HALFSPIN)
            : n(chk_principal_num(n))
            , l(chk_azimuthal_num(n, l))
            , m(chk_magnetic_num(l, m))
            , s(chk_spin_num(s))
        {
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Static methods
        ///

        static bool is_principal_num_ok(int n) { return 1 <= n; }
        static bool is_azimuthal_num_ok(int n, int l) { return 0 <= l && l <= n - 1; }
        static bool is_magnetic_num_ok(int l, int m) { return -l <= m && m <= l; }
        static bool is_spin_num_ok(double s) { return s == HALFSPIN || s == -HALFSPIN; }

        static int chk_principal_num(int n) { return is_principal_num_ok(n) ? n : throw invalid_argument("Invalid principal number n = " + to_string(n)); }
        static int chk_azimuthal_num(int n, int l) { return is_azimuthal_num_ok(n, l) ? l : throw invalid_argument("Invalid azimuthal number l = " + to_string(l)); }
        static int chk_magnetic_num(int l, int m) { return is_magnetic_num_ok(l, m) ? m : throw invalid_argument("Invalid magnetic number m = " + to_string(m)); }
        static double chk_spin_num(double s) { return is_spin_num_ok(s) ? s : throw invalid_argument("Invalid spin number s = " + to_string(s)); }
    };
}

#endif /* !QUANTUM_NUMBERS_HPP */