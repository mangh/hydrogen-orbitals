#ifndef MANGH_HYDROGEN_PARAMETERS_HPP
#define MANGH_HYDROGEN_PARAMETERS_HPP

#include "numbers.hpp"

using namespace mangh::math;

namespace mangh::hydrogen
{
    /**
     * @brief Physical parameters of the hydrogen atom.
     * 
     */
    struct Parameters
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Bohr radius.
        const number_t a0;

        /// @brief Planck constant.
        const number_t h;

        /// @brief Mass of particle.
        const number_t mp;

        /// @brief Potential constant (???)
        const number_t kp;

        /// @brief Constant in energy level equation (???).
        const number_t ke;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        Parameters() = delete;

        /**
         * @brief parameters constructor.
         * 
         * @param a0 - reduced (?) Bohr radius,
         * @param h  - Planck constant,
         * @param mp - Mass of particle,
         * @param kp - Potential constant (depending on the charge of the electron - ???),
         * @param ke - Constant in energy level equation (???).
         *
         */
        Parameters(number_t a0, number_t h, number_t mp, number_t kp, number_t ke)
            : a0(a0)
            , h(h)
            , mp(mp)
            , kp(kp)
            , ke(ke)
        {
        }
    };
}

#endif /* !MANGH_HYDROGEN_PARAMETERS_HPP */