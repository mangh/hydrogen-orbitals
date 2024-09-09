#ifndef SPECIAL_POLYNOMIALS_HPP
#define SPECIAL_POLYNOMIALS_HPP

#include <cassert>
#include <cmath>
#include <type_traits>
#include <stdexcept>

namespace SpecialPolynomial
{
    constexpr double PI = 3.141592653589793238462643383279502884L;

    /**
     * @brief Generalized (associated) Laguerre polynomial.
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

    /**
     * @brief Renormalized associated Legendre polynomial ~Pm,l(x).
     * 
     * @details Computes renormalized (for better convergence) associated Legendre polynomial.
     * See: Numerical Recipes in C++, 6.7. Spherical Harmonics, equation (6.7.8).
     * 
     * @param l - the degree of the polynomial,
     * @param m - the order of the polynomial,
     * @param x - the argument.
     * @returns Value of the renormalized associated Legendre polynomial ~Pm,l(x).
     * 
     * @exception invalid_argument - when the arguments are outside their allowed range.
     * 
     * @pre 0 <= m <= l
     * @pre -1 <= x <= 1
     * 
     */
    template<typename T>
    T assoc_legendre_normalized(unsigned l, unsigned m, T x)
    {
        static_assert(std::is_floating_point<T>::value,
            "assoc_legendre_normalized() requires a floating-point type parameter x."
        );

        if (m > l || std::abs(x) > 1.0)
            throw std::invalid_argument("assoc_legendre_normalized: bad arguments.");

        // compute ~Pm,m
        T Pmm = 1.0;
        if (m > 0)
        {
            T omx2 = (1.0 - x) * (1.0 + x);
            T fact = 1.0;
            for (unsigned i = 1; i <= m; i++)
            {
                Pmm *= omx2 * fact / (fact + 1.0);
                fact += 2.0;
            }
        }
        
        Pmm = std::sqrt((2 * m + 1) * Pmm / (4.0 * PI));

        if (m & 1)
            Pmm = -Pmm;

        if (l == m)
            return Pmm;

        // compute ~Pm,m+1
        T oldfact = std::sqrt(2.0 * m + 3.0);
        T Pmmp1 = x * oldfact * Pmm;
        if (l == m+1)
            return Pmmp1;

        // compute ~Pm,l (m+1 < l)
        T Pll = 0.0;
        for (unsigned ll = m + 2; ll <= l; ll++)
        {
            T fact = std::sqrt((4.0 * ll * ll - 1.0) / (ll * ll - m * m));
            Pll = (x * Pmmp1 - Pmm / oldfact) * fact;
            oldfact = fact;
            Pmm = Pmmp1;
            Pmmp1 = Pll;
        }
        return Pll;
    }

    /**
     * @brief Associated Legendre polynomial Pm,l(x).
     * See: Numerical Recipes in C++, 6.7. Spherical Harmonics, equation (6.7.4).
     * 
     * @param l - the degree of the polynomial,
     * @param m - the order of the polynomial,
     * @param x - the argument.
     * @returns Value of the associated Legendre polynomial Pm,l(x).
     * @exception invalid_argument - when the arguments are outside their allowed range.
     * 
     * @pre 0 <= m <= l
     * @pre -1 <= x <= 1
     * 
     * @remarks The function may overflow for m >= 80, or even sooner if l >> m.
     * 
     */
    template<typename T>
    T assoc_legendre(unsigned l, unsigned m, T x)
    {
        if (m > l || std::abs(x) > 1.0)
            throw std::invalid_argument("assoc_legendre: bad arguments.");

        T prod = 1.0;
        for (unsigned j = l - m + 1; j <= l + m; j++)
        {
            prod *= j;
        } 
        return std::sqrt(4.0 * PI * prod / (2 * l + 1)) * assoc_legendre_normalized(l,m,x);
    }
}

#endif /* !SPECIAL_POLYNOMIALS_HPP */