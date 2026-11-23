#ifndef MANGH_MATH_FACTORIAL_HPP
#define MANGH_MATH_FACTORIAL_HPP

#include <cassert>
#include <cmath>
#include <type_traits>

namespace mangh { namespace math
{
    /**
     * @brief Factorial (that can return large values up to overflow that is not signalled!).
     * 
     * @tparam T - floating point type of the return value.
     * @param n - integer argument (n >= 0).
     * @returns n! = factorial(n).
     * 
     * @note The function can result in overflow (just like any factorial), but this is not checked.
     * 
     */
    template<typename T>
    inline T factorial(unsigned n)
    {
        static_assert(
            std::is_floating_point<T>::value,
            "T factorial(n): function template requires a floating-point type parameter T."
        );

        return std::tgamma(static_cast<T>(n) + T{1});
    }
}}

#endif /* !MANGH_MATH_FACTORIAL_HPP */
