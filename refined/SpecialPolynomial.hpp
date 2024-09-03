#ifndef SPECIAL_POLYNOMIALS_HPP
#define SPECIAL_POLYNOMIALS_HPP

#include <cassert>
#include <cmath>
#include <type_traits>

namespace SpecialPolynomial
{
    template<typename N>
    unsigned long factorial(N n)
    {
        static_assert(
            std::numeric_limits<N>::is_integer &&
            !std::numeric_limits<N>::is_signed,
            "factorial(n) requires a parameter of type unsigned integer."
        );
        return (n <= 0) ? 1 : n * factorial(n - 1);
    }

    /**
     * @brief Legendre polynomial Pn(x).
     * 
     * @tparam T - floating point type of the argument.
     * @param n - the degree of the polynomial,
     * @param x - the argument.
     * @returns T 
     *      
     * @remarks Recursive definition (https://mathworld.wolfram.com/LegendrePolynomial.html):
     *      
     *          1    n/2+1  (-1)^k∙(2l-2k)!∙x^(n-2k)
     * Pn(x) = --- ∙   Σ    ------------------------
     *         2^n    k=0      k!∙(n-k)!∙(n-2k)!
     *      
     */
    template<typename T>
    T Pn(unsigned n, T x)
    {
        static_assert(std::is_floating_point<T>::value, "Legendre polynomial Pn(n,x) requires a floating-point type parameter x.");
        T total = 0.0;
        for (unsigned k = 0; k < n / 2; ++k)
        {
            total += pow(-1, k)
                   * pow(x, n - 2 * k)
                   * factorial(2 * n - 2 * k)
                   / factorial(k)
                   / factorial(n - k)
                   / factorial(n - 2 * k);
        }
        return total / pow(2.0, n);
    }

    template<typename T>
    T dm_Pl_dxm(T x, int m, int l)
    {
        assert(0 <= m);
        T dx = 0.0001;
        return
            (m == 0) ? Pn(l, x) :
            (m == 1) ? (Pn(l, x + dx) - Pn(l, x)) / dx :
            (dm_Pl_dxm(x + dx, m - 1, l) - dm_Pl_dxm(x, m - 1, l)) / dx;
    }

    /**
     * @brief Associated Legendre polynomial Pnm(n,m,x).
     *
     * @tparam T - floating point type of the argument.
     * @param n - the degree of the polynomial,
     * @param m - the order of the polynomial,
     * @param x - the argument.
     * @returns T
     */
    template<typename T>
    T Pnm(unsigned n, int m, T x)
    {
        static_assert(std::is_floating_point<T>::value, "Associated Legendre polynomial Pnm(n,m,x) requires a floating-point type parameter x.");
        return(m >= 0) ?
            pow(-1, m) * pow(1.0 - x * x, m / 2.0) * dm_Pl_dxm(x, m, n) :
            pow(-1, (-m)) * factorial(n - (-m)) / factorial(n + (-m)) * Pnm(n, (-m), x);
    }

    /**
     * @brief Generalized (or associated) Laguerre polynomial Ln(n,a,x).
     * 
     * @tparam T - floating point type of the argument.
     * @param n - the degree of the polynomial,
     * @param a - the order of the polynomial,
     * @param x - the argument.
     * @returns T 
     * 
     * @remarks Recursive definition (https://en.wikipedia.org/wiki/Laguerre_polynomials):
     *      
     * L[0][a](x) = 1
     * L[1][a](x) = 1 + a -x
     *                (2k+1+a-x)∙L[n][a] - (k+a)∙L[n-1][a]
     * L[n+1][a](x) = ------------------------------------
     *                                n + 1
     */
    template<typename T>
    T Ln(unsigned n, int a, T x)
    {
        static_assert(std::is_floating_point<T>::value, "Associated Laguerre polynomial Ln(n,a,x) requires a floating-point type parameter x.");
        return
            (n == 0) ? 1.0 :
            (n == 1) ? 1.0 + a - x :
            ((2.0 * (n - 1) + 1.0 + a - x) * Ln(n - 1, a, x) - (n - 1 + a) * Ln(n - 2, a, x)) / n;
    }
}

#endif /* !SPECIAL_POLYNOMIALS_HPP */