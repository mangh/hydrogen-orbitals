#ifndef SPECIAL_POLYNOMIALS_HPP
#define SPECIAL_POLYNOMIALS_HPP

#include <cassert>
#include <cmath>
#include <type_traits>

namespace SpecialPolynomial
{
    /**
     * @brief Factorial that can return large values.
     * 
     * @tparam T - floating point type of the argument.
     * @param n - factorial argument.
     * @returns Factorial value i.e. n!
     */
    template<typename T>
    T factorial(T n)
    {
        static_assert(std::is_floating_point<T>::value, "factorial(n) requires a floating-point type parameter.");
        return std::exp(std::lgamma(n + 1));
    }

    /**
     * @brief Binomial coefficient "n choose k" (n/k).
     * 
     * @param n - numerator argument,
     * @param k - denominator argument.
     * @returns The coefficient (n/k) in the expansion of the polynomial (x+y)^n.
     * 
     * @remarks https://stackoverflow.com/questions/44718971/calculate-binomial-coffeficient-very-reliably
     */
    template<typename T>
    T binom(T n, T k)
    {
        static_assert(std::is_floating_point<T>::value, "binom(n, k) requires a floating-point type parameters.");
        return std::exp(
            std::lgamma(n + 1) - std::lgamma(n - k + 1) - std::lgamma(k + 1)
        );
    }

    /**
     * @brief Legendre polynomial.
     * 
     * @tparam T - floating point type of the argument.
     * @param n - the degree of the polynomial,
     * @param x - the argument.
     * @returns T 
     *      
     * @remarks Recursive definition (https://mathworld.wolfram.com/LegendrePolynomial.html):
     *      
     *                1    n/2+1           /n\     /2n-2k\
     * legendre(x) = --- *   Σ   (-1)^k * |   | * |       | * x^(n-2k)
     *               2^n    k=0            \k/     \  n  /
     *      
     *                1    n/2+1              (2n-2k)!
     *             = --- *   Σ   (-1)^k * ----------------- * x^(n-2k)
     *               2^n    k=0           k!*(n-k)!*(n-2k)!
     *      
     */
    template<typename T>
    T legendre(unsigned n, T x)
    {
        static_assert(std::is_floating_point<T>::value, "Legendre polynomial requires a floating-point type parameter x.");
        T total = 0.0L;
        for (unsigned k = 0; k < n / 2; ++k)
        {
            // total += pow(-1, k)
            //        * pow(x, n - 2 * k)
            //        * factorial(2 * n - 2 * k)
            //        / factorial(k)
            //        / factorial(n - k)
            //        / factorial(n - 2 * k);
            total += pow(-1, k)
                   * pow(x, n - 2 * k)
                   * binom(static_cast<double>(n), static_cast<double>(k))
                   * binom(static_cast<double>(2 * n - 2 * k), static_cast<double>(n));
        }
        return total / pow(2.0, n);
    }

    template<typename T>
    T dm_Pl_dxm(T x, int m, int l)
    {
        assert(0 <= m);
        T dx = 0.0001L;
        return
            (m == 0) ? legendre(l, x) :
            (m == 1) ? (legendre(l, x + dx) - legendre(l, x)) / dx :
            (dm_Pl_dxm(x + dx, m - 1, l) - dm_Pl_dxm(x, m - 1, l)) / dx;
    }

    /**
     * @brief Associated Legendre polynomial.
     *
     * @tparam T - floating point type of the argument.
     * @param n - the degree of the polynomial,
     * @param m - the order of the polynomial,
     * @param x - the argument.
     * @returns T
     */
    template<typename T>
    T assoc_legendre(unsigned n, int m, T x)
    {
        static_assert(std::is_floating_point<T>::value, "Associated Legendre polynomial assoc_legendre(n,m,x) requires a floating-point type parameter x.");
         return(m >= 0) ?
             pow(-1, m) * pow(1.0 - x * x, m / 2.0) * dm_Pl_dxm(x, m, n) :
             pow(-1, (-m))
                * factorial(static_cast<double>(n - (-m)))
                / factorial(static_cast<double>(n + (-m)))
                * assoc_legendre(n, (-m), x);
    }

    /**
     * @brief Generalized (or associated) Laguerre polynomial.
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
    T assoc_laguerre(unsigned n, int a, T x)
    {
        static_assert(std::is_floating_point<T>::value, "Associated Laguerre polynomial requires a floating-point type parameter x.");
        return
            (n == 0) ? 1.0 :
            (n == 1) ? 1.0 + a - x :
            ((2.0 * (n - 1) + 1.0 + a - x)
                * assoc_laguerre(n - 1, a, x) - (n - 1 + a)
                * assoc_laguerre(n - 2, a, x)) / n;
    }
}

#endif /* !SPECIAL_POLYNOMIALS_HPP */