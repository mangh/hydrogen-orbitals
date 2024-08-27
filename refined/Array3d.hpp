#ifndef ARRAY3D_HPP
#define ARRAY3D_HPP

#include <cassert>
#include <complex>
#include <vector>

#include "cppfeatures.h"

#include "DomainData.hpp"

namespace HydrogenOrbitals
{
    using ComplexData = std::complex<double>;

    /**
     * @brief 3-dimensional array (A[i][j][k]) of numbers (including complex numbers)
     * 
     */
    template<typename T = std::complex<double>>
    struct Array3d
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Numeric (complex) data
        std::vector<
            std::vector<
                std::vector<T>>> m_arry;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        Array3d() = delete;

        /**
         * @brief Array3d constructor.
         * 
         * @param size_i - the number of planes,
         * @param size_j - the number of rows in a single plane,
         * @param size_k - the number of columns in a single row.
         * 
         */
        Array3d(std::size_t size_i, std::size_t size_j, std::size_t size_k)
            : m_arry(size_i)
        {
            for(std::size_t i = 0; i < size_i; ++i)
            {
                m_arry[i] = std::vector<std::vector<T>>(size_j);
                for(std::size_t j = 0; j < size_j; ++j)
                {
                    m_arry[i][j] = std::vector<T>(size_k);
                }
            }
        }

        /**
         * @brief Array3d constructor.
         * @param dd - Domain Data.
         */
        Array3d(const DomainData& dd)
            : Array3d(dd.n_r, dd.n_theta, dd.n_phi)
        {
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Operators
        ///

        /**
         * @brief Array subscript operator
         * 
         * @param r 
         * @returns constexpr std::vector<std::vector<ComplexData>>& 
         */
        CONSTEXPR_SINCE_CXX20 std::vector<std::vector<T>>& operator[](std::size_t i) const
        {
            assert(i < m_arry.size());
            return m_arry[i];
        }
        std::vector<std::vector<T>>& operator[](std::size_t i)
        {
            assert(i < m_arry.size());
            return m_arry[i];
        }

        CONSTEXPR_SINCE_CXX20 T& operator()(std::size_t i, std::size_t j, std::size_t k) const
        {
            assert(i < m_arry.size() && j < m_arry[i].size() && k < m_arry[i][j].size());
            return m_arry[i][j][k];
        }

        T& operator()(std::size_t i, std::size_t j, std::size_t k)
        {
            assert(i < m_arry.size() && j < m_arry[i].size() && k < m_arry[i][j].size());
            return m_arry[i][j][k];
        }
    };
}

#endif /* !ARRAY3D_HPP */