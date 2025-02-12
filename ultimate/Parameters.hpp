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

        /// @brief Planck constant.
        const number_t h;

        /// @brief Mass of particle.
        const number_t mp;

        /// @brief Potential energy constant.
        ///
        /// Term for the potential energy of the electron in the electrostatic field of the proton:
        ///
        ///   - V(r) = -1 / (4·π·ε0) · Z·e^2 / r
        ///
        /// where:
        /// 
        ///   - ε0 is the permittivity of the vacuum,
        /// 
        ///   - Z is the atomic number (number of protons in the nucleus),
        /// 
        ///   - e is the elementary charge (charge of an electron),
        /// 
        ///   - r is the distance of the electron from the nucleus.
        /// 
        const number_t kp;

        /// @brief Energy constant (negative value of Rydberg's energy unit).
        ///
        /// The energy levels of hydrogen-like atoms:
        ///
        ///   - En = -Rh / n^2
        ///
        /// where:
        ///
        ///   - n is the principal quantum number,
        ///
        ///   - Rh is the Rydberg's energy unit (Rh = 2.18·10^-18 J).
        ///
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
        Parameters(number_t h, number_t mp, number_t kp, number_t ke)
            : h(h)
            , mp(mp)
            , kp(kp)
            , ke(ke)
        {
        }
    };
}

#endif /* !MANGH_HYDROGEN_PARAMETERS_HPP */