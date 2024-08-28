#ifndef PHYSICAL_PARAMETERS_HPP
#define PHYSICAL_PARAMETERS_HPP

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
        const double a0;

        /// @brief Planck (???) constant.
        const double h;

        /// @brief Mass of particle.
        const double mp;

        /// @brief Potential constant (???)
        const double kp;

        /// @brief Constant in energy level equation (???).
        const double ke;

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
        PhysicalParameters(double a0, double h, double mp, double kp, double ke)
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