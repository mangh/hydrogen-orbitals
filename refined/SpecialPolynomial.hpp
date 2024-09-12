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
    * @brief Factorial that can return large values.
    * 
    * @tparam T - floating point type of the return value.
    * @param n - factorial argument.
    * @returns Factorial value i.e. n!
    */
    template<typename T>
    static T fact(unsigned n)
    {
        static_assert(std::is_floating_point<T>::value, "fact(n) requires a floating-point type parameter.");
        return std::exp(std::lgamma(n + 1));
    }

    /**
     * @brief Generalized (associated) Laguerre polynomial.
     * 
     * @tparam T - floating point type of the argument.
     * 
     * @param n - the degree of the polynomial,
     * @param alpha - the order of the polynomial,
     * @param x - the argument.
     * @returns T 
     *
     * @remarks See: Numerical Recipes in C++, 4.6.1, algorithm "gauleg",
     * and the recursive definition at https://en.wikipedia.org/wiki/Laguerre_polynomials:
     * 
     * L[alpha][0](x) = 1
     * L[alpha][1](x) = 1 + alpha - x
     * L[alpha][k+1](x) = ((2k+1+alpha-x) * L[alpha][k](x) - (k+alpha) * L[alpha][k-1](x)) / (k + 1).
     * 
     */
    template<typename T>
    T assoc_laguerre(unsigned n, unsigned alpha, T x)
    {
        static_assert(std::is_floating_point<T>::value,
            "Associated Laguerre polynomial requires a floating-point type parameter x."
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

    /**
     * @brief Associated Legendre polynomial normalized for use in spherical harmonics.
     * 
     * @tparam T - floating point type of the argument.
     * 
     * @param l - the degree of the polynomial,
     * @param m - the order of the polynomial,
     * @param x - the argument.
     * 
     * @returns The value of the associated Legendre polynomial ~Pm,l(x), however, normalized for use with spherical harmonics Yl,m.
     * 
     * @pre 0 <= m <= l
     * @pre -1 <= x <= 1
     * @exception invalid_argument - when the arguments are outside their allowed range.
     * 
     * @remarks
     * See: Numerical Recipes in C++, 6.7. Spherical Harmonics, equation (6.7.8), algorithm "plegendre".
     * It converges better (than the algorithm for the standard Legendre polynomial).
     */
    template<typename T>
    T ylm_assoc_legendre(unsigned l, unsigned m, T x)
    {
        static_assert(std::is_floating_point<T>::value,
            "ylm_assoc_legendre() requires a floating-point type parameter x."
        );

        if (m > l)
            throw std::invalid_argument("ylm_assoc_legendre(l,m,x): bad l, m arguments. Should be: m <= l, and is: m > l.");
        if ((x < -1) || (x > 1))
            throw std::invalid_argument("ylm_assoc_legendre(l,m,x): bad argument x. Should be: -1 <= x <= 1, but it is not.");

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
     * @brief Factor to renormalize Legendre polynomials built for spherical harmonics.
     * 
     * @tparam T - floating point type of the argument.
     * 
     * @param l - the degree of the polynomial,
     * @param m - the order of the polynomial,
     * 
     * @returns The factor: 
     * 
     *    sqrt(((2l+1) / (4pi)) * ((l-m)! / (l+m)!))
     * 
     * required to renormalize the ylm_assoc_legendre() polynomial
     * to the standard (orthonormal) Legendre polynomial:
     * 
     *    assoc_legendre(l,m,x) = 
     *        ylm_assoc_legendre(l,m,x) / ylm_normalizing_factor<T>(l,m).
     * 
     * @pre m <= l
     * @exception invalid_argument - when the @precondition is not met.
     * 
     */
    template<typename T>
    T ylm_normalizing_factor(unsigned l, unsigned m)
    {
        static_assert(std::is_floating_point<T>::value,
            "ylm_assoc_legendre() requires a floating-point type parameter x."
        );

        if (m > l)
            throw std::invalid_argument("ylm_normalizing_factor(l,m): bad arguments. Should be: m <= l, and is: m > l.");

        T prod = 1.0;
        for (unsigned j = l - m + 1; j <= l + m; ++j) prod *= j;
        return std::sqrt((2.0 * l + 1) / 4.0 / PI / prod);
    }

    /**
     * @brief Associated Legendre polynomial Pm,l(x) (orthonormal).
     * 
     * @tparam T - floating point type of the argument.
     * 
     * @param l - the degree of the polynomial,
     * @param m - the order of the polynomial,
     * @param x - the argument.
     * @returns Value of the associated Legendre polynomial Pm,l(x).
     * 
     * @pre 0 <= m <= l
     * @pre -1 <= x <= 1
     * @exception invalid_argument - when the arguments are outside their allowed range.
     * 
     * @remarks See: Numerical Recipes in C++, 6.7. Spherical Harmonics, equation (6.7.4), algorithm "plgndr".
     * The function may overflow for m >= 80, or even sooner if l >> m.
     * 
     */
    template<typename T>
    T assoc_legendre(unsigned l, unsigned m, T x)
    {
        static_assert(std::is_floating_point<T>::value,
            "ylm_assoc_legendre() requires a floating-point type parameter x."
        );
        return ylm_assoc_legendre(l,m,x) / ylm_normalizing_factor<T>(l, m);
    }
}

#endif /* !SPECIAL_POLYNOMIALS_HPP */