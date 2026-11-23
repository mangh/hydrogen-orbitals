#ifndef MANGH_LINSPACE_HPP
#define MANGH_LINSPACE_HPP

#include <cstddef>
#include <type_traits>
#include <vector>

namespace mangh { namespace math
{
    /**
    * @brief linspace (a vector of evenly spaced numbers of type T).
    *
    */
    template<typename T>
    struct linspace : std::vector<T>
    {
        static_assert(std::is_floating_point_v<T>, "linspace<T>: T must be a floating-point type.");

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        /// @brief Delta (step) value - the interval between consecutive numbers.
        const T delta;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /// @brief No default constructor!
        linspace() = delete;

        /**
        * @brief linspace constructor.
        *
        * @param start - the start value,
        * @param end   - the end value,
        * @param size  - the number of points between the `start` and `end` values.
        *
        * @retval
        *
        * @li - for `size > 1`: vector of `size` elements evenly spaced between the `start` and `end` values
        *       (the last element is exactly `end`),
        *
        * @li - for `size = 1`: one-element vector (with the only element set to `start`),
        *
        * @li - for `size < 1`: empty vector.
        *
        */
        linspace(T start, T end, std::size_t size)
            : delta(size > 1 ? (end - start) / static_cast<T>(size - 1) : T{})
        {
            this->reserve(size);
            for (std::size_t i = 0; i < size; ++i) this->push_back(start + static_cast<T>(i) * delta);
            if (size > 1) this->back() = end;
        }
    };
}}

#endif /* !MANGH_LINSPACE_HPP */
