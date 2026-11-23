#ifndef MANGH_ARRAY3D_HPP
#define MANGH_ARRAY3D_HPP

#include <cassert>
#include <cstddef>
#include <vector>

#include "cppfeatures.hpp"

namespace mangh { namespace math
{
    /**
     * @brief 3-dimensional array stored contiguously in C (row-major) order.
     *
     * The element (p, r, c) is at the offset `(p · rows + r) · columns + c`,
     * so the whole array can be written at once (e.g. to a .npy file, see npy::save).
     *
     */
    template<typename T>
    class array3d
    {
    public:

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        array3d() = delete;

        /**
         * @brief array3d constructor (all elements value-initialized).
         *
         * @param pages   - the number of pages,
         * @param rows    - the number of rows (in a single page),
         * @param columns - the number of columns (in a single row).
         *
         */
        array3d(std::size_t pages, std::size_t rows, std::size_t columns)
            : m_pages(pages)
            , m_rows(rows)
            , m_columns(columns)
            , m_data(pages * rows * columns)
        {
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Operators
        ///

        /**
         * @brief Function-call operator to access a single data item (for reading).
         *
         * @param p - page index,
         * @param r - row index,
         * @param c - column index.
         * @returns Single data item (read-only) for the given indexes (p, r, c).
         */
        CONSTEXPR_SINCE_CXX20 const T& operator()(std::size_t p, std::size_t r, std::size_t c) const
        {
            return m_data[offset(p, r, c)];
        }

        /**
         * @brief Function-call operator to access a single data item (for writing).
         *
         * @param p - page index,
         * @param r - row index,
         * @param c - column index.
         * @returns Single data item (writable) for the given indexes (p, r, c).
         */
        CONSTEXPR_SINCE_CXX20 T& operator()(std::size_t p, std::size_t r, std::size_t c)
        {
            return m_data[offset(p, r, c)];
        }

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Methods
        ///
        ///     (std::vector is usable in constant expressions since C++20, hence CONSTEXPR_SINCE_CXX20)
        ///

        constexpr std::size_t pages() const noexcept { return m_pages; }
        constexpr std::size_t rows() const noexcept { return m_rows; }
        constexpr std::size_t columns() const noexcept { return m_columns; }

        /// @brief The total number of elements (pages · rows · columns).
        CONSTEXPR_SINCE_CXX20 std::size_t size() const noexcept { return m_data.size(); }

        /// @brief Pointer to the first element (the elements are contiguous, in C order).
        CONSTEXPR_SINCE_CXX20 const T* data() const noexcept { return m_data.data(); }
        CONSTEXPR_SINCE_CXX20 T* data() noexcept { return m_data.data(); }

    private:

        CONSTEXPR_SINCE_CXX20 std::size_t offset(std::size_t p, std::size_t r, std::size_t c) const
        {
            assert(p < m_pages && r < m_rows && c < m_columns);
            return (p * m_rows + r) * m_columns + c;
        }

        std::size_t m_pages;
        std::size_t m_rows;
        std::size_t m_columns;

        /// @brief Data pages (each page consists of rows, and each row consists of columns), one after another.
        std::vector<T> m_data;
    };
}}

#endif /* !MANGH_ARRAY3D_HPP */
