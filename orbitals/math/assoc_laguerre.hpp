#ifndef MANGH_MATH_ASSOC_LAGUERRE_HPP
#define MANGH_MATH_ASSOC_LAGUERRE_HPP

#include <cassert>
#include <type_traits>

namespace mangh { namespace math
{
    /**
     * @brief Associated (generalized) Laguerre polynomial `Lα,n(x)`.
     * 
     * @tparam T - floating point type of the argument.
     * 
     * @param alpha - the parameter α of the polynomial (α >= 0),
     * @param n - the degree of the polynomial (n >= 0),
     * @param x - the argument of the polynomial.
     * 
     * @returns The value of the polynomial `Lα,n(x)` at `x`.
     *
     * @remarks See: Numerical Recipes in C++, 4.6.1, algorithm "gauleg",
     * and the recursive definition at https://en.wikipedia.org/wiki/Laguerre_polynomials:
     * 
     * Lα,0(x) = 1,
     * Lα,1(x) = 1 + α - x,
     * Lα,k+1(x) = ((2k+1+α-x) · Lα,k(x) - (k+α) · Lα,k-1(x)) / (k + 1).
     * 
     */
    template<typename T>
    T assoc_laguerre(unsigned alpha, unsigned n, T x)
    {
        static_assert(
            std::is_floating_point<T>::value,
            "T assoc_laguerre(unsigned n, unsigned alpha, T x): template requires a floating-point type parameter T."
        );

        T L2 = 0.0;
        T L1 = 1.0;
        for(unsigned i = 0; i < n; ++i)
        {
            T L3 = L2;
            L2 = L1;
            L1 = ((2 * i + 1 + alpha - x) * L2 - (i + alpha) * L3) / (i + 1);
        }
        return L1;
    }
}}

#endif /* !MANGH_MATH_ASSOC_LAGUERRE_HPP */
