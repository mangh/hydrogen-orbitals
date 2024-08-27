#ifndef DOMAIN_DATA_HPP
#define DOMAIN_DATA_HPP

#include <cstddef>

namespace HydrogenOrbitals
{
    /**
    * @brief Parameters
    * 
    */
    struct DomainData
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Number of nodes in r direction.
        const std::size_t n_r;

        /// @brief Number of nodes in theta direction.
        const std::size_t n_theta;

        /// @brief Number of nodes in phi direction.
        const std::size_t n_phi;

        /// @brief Length of domain.
        const double R;

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
        * @param ct_r     - number of nodes in r direction,
        * @param ct_theta - number of nodes in theta direction
        * @param ct_phi   - number of nodes in phi direction
        * @param R_length - length of domain
        *
        * @remarks Data used in Derek's (original) version: n_r = 200, n_theta = 40, n_phi = 20, R = 40.0.
        *
        */
        DomainData(std::size_t ct_r, std::size_t ct_theta, std::size_t ct_phi, double R_length)
            : n_r(ct_r)
            , n_theta(ct_theta)
            , n_phi(ct_phi)
            , R(R_length)
        {
        }
    };
}

#endif /* !DOMAIN_DATA_HPP */