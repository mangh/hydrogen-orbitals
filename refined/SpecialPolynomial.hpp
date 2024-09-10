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
     * 
     * L[1][a](x) = 1 + a - x
     * 
     * L[k+1][a](x) = ((2k+1+a-x) * L[k][a](x) - (k+a) * L[k-1][a](x)) / (k + 1).
     */
    template<typename T>
    T assoc_laguerre(unsigned n, unsigned a, T x)
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
     * @brief Associated Legendre polynomial for use with spherical harmonics.
     * 
     * @details Computes the associated Legendre polynomial ~Pm,l(x) prenormalized for orthonormal spherical harmonics Yl,m.
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
     * @remarks See: Numerical Recipes in C++, 6.7. Spherical Harmonics, equation (6.7.8), algorithm "plegendre".
     * The square root factor used in equation (6.7.2) makes the resulting polynomial ~Pl,m well suited for spherical harmonics,
     * but, importantly, it also prevents the algorithm from diverging too early.
     * 
     */
    template<typename T>
    T assoc_legendre_ylm(unsigned l, unsigned m, T x)
    {
        static_assert(std::is_floating_point<T>::value,
            "assoc_legendre_ylm() requires a floating-point type parameter x."
        );

        if (m > l)
            throw std::invalid_argument("assoc_legendre_ylm(l,m,x): bad l, m arguments. Should be: m <= l, and is: m > l.");
        if ((x < -1) || (x > 1))
            throw std::invalid_argument("assoc_legendre_ylm(l,m,x): bad argument x. Should be: -1 <= x <= 1, but it is not.");

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
     * @remarks See: Numerical Recipes in C++, 6.7. Spherical Harmonics, equation (6.7.4), algorithm "plgndr".
     * The function may overflow for m >= 80, or even sooner if l >> m.
     * 
     */
    template<typename T>
    T assoc_legendre(unsigned l, unsigned m, T x)
    {
        if (m > l)
            throw std::invalid_argument("assoc_legendre(l,m,x): bad l, m arguments. Should be: m <= l, and is: m > l.");
        if ((x < -1) || (x > 1))
            throw std::invalid_argument("assoc_legendre(l,m,x): bad argument x. Should be: -1 <= x <= 1, but it is not.");

        T prod = 1.0;
        for (unsigned j = l - m + 1; j <= l + m; j++)
        {
            prod *= j;
        } 
        return std::sqrt(4.0 * PI * prod / (2 * l + 1)) * assoc_legendre_ylm(l,m,x);
    }
}

#endif /* !SPECIAL_POLYNOMIALS_HPP */