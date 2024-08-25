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
         * @param a0_ - Bohr radius,
         * @param h_  - Planck (???) constant,
         * @param mp_ - Mass of particle,
         * @param kp_ - Potential constant (???),
         * @param ke_ - Constant in energy level equation (???).
         *
         * @remarks Data used in Derek's (original) version: a0 = 1.0, h = 1.0, mp = 1.0, kp = 1.0, ke = 1.0.
         *
         */
        PhysicalParameters(double a0_, double h_, double mp_, double kp_, double ke_)
            : a0(a0_)
            , h(h_)
            , mp(mp_)
            , kp(kp_)
            , ke(ke_)
        {
        }
    };
}

#endif /* !PHYSICAL_PARAMETERS_HPP */