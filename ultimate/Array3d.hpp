#ifndef MANGH_ARRAY3D_HPP
#define MANGH_ARRAY3D_HPP

#include <cassert>
#include <vector>

#include "cppfeatures.hpp"

namespace mangh { namespace math
{
    /**
     * @brief 3-dimensional array.
     * 
     */
    template<typename T>
    struct Array3d
    {
        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Data pages(each page consists of rows, and each row consists of columns).
        std::vector<
            std::vector<
                std::vector<T>>> m_page;

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
         * @param pages  - the number of pages,
         * @param rows    - the number of rows (in a single page),
         * @param columns - the number of columns (in a single row).
         * 
         */
        Array3d(std::size_t pages, std::size_t rows, std::size_t columns)
            : m_page(pages)
        {
            for(std::size_t p = 0; p < pages; ++p)
            {
                m_page[p] = std::vector<std::vector<T>>(rows);
                for(std::size_t r = 0; r < rows; ++r)
                {
                    m_page[p][r] = std::vector<T>(columns);
                }
            }
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Operators
        ///

        /**
         * @brief Array subscript operator (for reading).
         * 
         * @param p - page index. 
         * @returns Data page (read-only) for the p-index.
         */
        CONSTEXPR_SINCE_CXX20 std::vector<std::vector<T>>& operator[](std::size_t p) const
        {
            assert(p < m_page.size());
            return m_page[p];
        }
        /**
         * @brief Array subscript operator (for writing).
         * 
         * @param p - page index. 
         * @returns Data page (for writing) for the p-index.
         */
        std::vector<std::vector<T>>& operator[](std::size_t p)
        {
            assert(p < m_page.size());
            return m_page[p];
        }

        /**
         * @brief Function-call operator to acces single data item (for reading).
         * 
         * @param p - page index,
         * @param r - row index,
         * @param c - column index.
         * @returns Single data item (read-only) for the given indexes (p, r, c).
         */
        CONSTEXPR_SINCE_CXX20 T& operator()(std::size_t p, std::size_t r, std::size_t c) const
        {
            assert(p < m_page.size() && r < m_page[p].size() && c < m_page[p][r].size());
            return m_page[p][r][c];
        }

        /**
         * @brief Function-call operator to acces single data item (for writing).
         * 
         * @param p - page index,
         * @param r - row index,
         * @param c - column index.
         * @returns Single data item (writable) for the given indexes (p, r, c).
         */
        T& operator()(std::size_t p, std::size_t r, std::size_t c)
        {
            assert(p < m_page.size() && r < m_page[p].size() && c < m_page[p][r].size());
            return m_page[p][r][c];
        }
    };
}}

#endif /* !MANGH_ARRAY3D_HPP */