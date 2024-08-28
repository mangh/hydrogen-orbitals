#ifndef PHYSICAL_PARAMETERS_HPP
#define PHYSICAL_PARAMETERS_HPP

#include "Numeral.hpp"
namespace HydrogenOrbitals
{
    /**
     * @brief Physical parameters.
     * 
     */
    struct PhysicalParameters
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Bohr radius.
        const numeral_t a0;

        /// @brief Planck (???) constant.
        const numeral_t h;

        /// @brief Mass of particle.
        const numeral_t mp;

        /// @brief Potential constant (???)
        const numeral_t kp;

        /// @brief Constant in energy level equation (???).
        const numeral_t ke;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        PhysicalParameters() = delete;

        /**
         * @brief PhysicalParameters constructor.
         * 
         * @param a0 - Bohr radius,
         * @param h  - Planck (???) constant,
         * @param mp - Mass of particle,
         * @param kp - Potential constant (depending on the charge of the electron - ???),
         * @param ke - Constant in energy level equation (???).
         *
         * @remarks Data used in Derek's (original) version: a0 = 1.0, h = 1.0, mp = 1.0, kp = 1.0, ke = 1.0.
         *
         */
        PhysicalParameters(numeral_t a0, numeral_t h, numeral_t mp, numeral_t kp, numeral_t ke)
            : a0(a0)
            , h(h)
            , mp(mp)
            , kp(kp)
            , ke(ke)
        {
        }
    };
}

#endif /* !PHYSICAL_PARAMETERS_HPP */