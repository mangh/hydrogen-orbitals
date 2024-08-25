#ifndef QUANTUM_NUMBERS_HPP
#define QUANTUM_NUMBERS_HPP

namespace HydrogenOrbitals
{
    /**
     * @brief Quantum numbers.
     * 
     */
    struct QuantumNumbers
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Principal quantum number.
        const int n;

        /// @brief Azimuthal number.
        const int l;

        /// @brief Magnetic quantum number.
        const int m;

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
         * @param n_ - principal quantum number,
         * @param l_ - azimuthal number,
         * @param m_ - magnetic quantum number.
         *
         * @remarks Data used in Derek's (original) version: n = 4, l = 2, m = 1.
         *
         */
        QuantumNumbers(int n_, int l_, int m_)
            : n(n_)
            , l(l_)
            , m(m_)
        {
        }
    };
}

#endif /* !QUANTUM_NUMBERS_HPP */