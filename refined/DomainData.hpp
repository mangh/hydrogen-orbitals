#ifndef DOMAIN_DATA_HPP
#define DOMAIN_DATA_HPP

#include <cstddef>
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

        /// @brief Number of nodes in r-direction (radial coordinate).
        const std::size_t n_r;

        /// @brief Number of nodes in θ-direction (azimuthal angle relative to the X-axis).
        const std::size_t n_theta;

        /// @brief Number of nodes in φ-direction (polar angle relative to the vertical Z-axis).
        const std::size_t n_phi;

        /// @brief Length of domain.
        const numeral_t R;

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
        * @remarks Data used in Derek's (original) version: n_r = 200, n_theta = 40, n_phi = 20, R = 40.0.
        *
        */
        DomainData(std::size_t n_r, std::size_t n_theta, std::size_t n_phi, numeral_t R)
            : n_r(n_r)
            , n_theta(n_theta)
            , n_phi(n_phi)
            , R(R)
        {
        }
    };
}

#endif /* !DOMAIN_DATA_HPP */