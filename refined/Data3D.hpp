#ifndef MAT3D_HPP
#define MAT3D_HPP

#include <cassert>
#include <complex>
#include <vector>
#include "DomainData.hpp"

namespace HydrogenOrbitals
{
    using ComplexData = std::complex<double>;

    /**
     * @brief Data3D
     * 
     */
    struct Data3D
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Numeric (complex) data
        std::vector<
            std::vector<
                std::vector<ComplexData>>> m_data;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        Data3D() = delete;

        /**
         * @brief Data3D constructor.
         * @param dd - Domain Data.
         */
        Data3D(const DomainData& dd)
            : m_data(dd.n_r)
        {
            for(int i = 0; i < dd.n_r; ++i)
            {
                m_data[i] = std::vector<std::vector<ComplexData>>(dd.n_theta);
                for(int j = 0; j < dd.n_theta; ++j)
                {
                    m_data[i][j] = std::vector<ComplexData>(dd.n_phi);
                }
            }
            // for(auto &r_dat : m_data)
            // {
            //     r_dat = std::vector<std::vector<ComplexData>>(dd.n_theta);
            //     for(auto &theta_dat : r_dat)
            //     {
            //         theta_dat = std::vector<ComplexData>(dd.n_phi);
            //     }
            // }
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Operators
        ///

        // constexpr ComplexFloat& operator[](std::size_t z, std::size_t y, std::size_t x) // C++23
        // {
        //     assert(x < X and y < Y and z < Z);
        //     return m[z * Y * X + y * X + x];
        // }

        /**
         * @brief Array subscript operator
         * 
         * @param r 
         * @returns constexpr std::vector<std::vector<ComplexData>>& 
         */
        constexpr std::vector<std::vector<ComplexData>>& operator[](std::size_t r)
        {
            assert(r < m_data.size());
            return m_data[r];
        }

        constexpr ComplexData& operator()(std::size_t r, std::size_t theta, std::size_t phi)
        {
            assert(r < m_data.size() && theta < m_data[r].size() && phi < m_data[r][theta].size());
            return m_data[r][theta][phi];
        }
    };
}

#endif /* !MAT3D_HPP */