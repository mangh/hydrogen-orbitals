#ifndef DOMAIN_DATA_HPP
#define DOMAIN_DATA_HPP

#include <cassert>
#include "Numeral.hpp"

namespace HydrogenOrbitals
{
    /**
    * @brief 3D grid parameters.
    * 
    */
    struct DomainData
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///
        ///     Notes: spherical coordinates (r, θ, φ)
        ///     by mathematical (not physical) convention:
        ///
        ///       r ≥ 0: radial coordinate,
        ///       0° ≤ θ < 360°: azimuthal angle (in the XY plane, measured from the X axis in the counterclockwise direction),
        ///       0° ≤ φ ≤ 180°: polar angle (the reference angle from z-axis).
        ///
        ///     See: https://en.wikipedia.org/wiki/Spherical_coordinate_system
        ///

        /// @brief Length of domain.
        const numeral_t R;

        /// @brief Number of nodes in r-direction (radial coordinate).
        const std::size_t n_r;

        /// @brief Number of nodes in θ-direction
        const std::size_t n_theta;

        /// @brief Number of nodes in φ-direction.
        const std::size_t n_phi;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        DomainData() = delete;

        /**
        * @brief DomainData constructor
        * 
        * @param n_r     - number of nodes in "r" direction,
        * @param n_theta - number of nodes in "theta" direction,
        * @param n_phi   - number of nodes in "phi" direction,
        * @param R       - length of domain.
        *
        */
        DomainData(numeral_t R, std::size_t n_r, std::size_t n_theta, std::size_t n_phi)
            : R(R)
            , n_r(n_r)
            , n_theta(n_theta)
            , n_phi(n_phi)
        {
            assert(R > 0.0);
            assert(n_r > 0);
            assert(n_theta > 0);
            assert(n_phi > 0);
        }
    };
}

#endif /* !DOMAIN_DATA_HPP */