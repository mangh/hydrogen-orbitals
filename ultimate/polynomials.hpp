#ifndef MANGH_POLYNOMIALS_HPP
#define MANGH_POLYNOMIALS_HPP

#include <cassert>
#include <cmath>
#include <type_traits>
#include <stdexcept>

#include "numbers.hpp"

namespace mangh::math
{
    /**
    * @brief Factorial (that can return large values).
    * 
    * @tparam T - floating point type of the return value.
    * @param n - integer argument (n >= 0).
    * @returns n! = factorial(n).
    */
    template<typename T>
    static T factorial(unsigned n)
    {
        static_assert(
            std::is_floating_point<T>::value,
            "T factorial(n): function template requires a floating-point type parameter T."
        );
        return std::exp(std::lgamma(n + 1));
    }

    /**
     * @brief Associated (generalized) Laguerre polynomial.
     * 
     * @tparam T - floating point type of the argument.
     * 
     * @param alpha - the parameter α of the polynomial (α >= 0),
     * @param n - the degree of the polynomial (n >= 0),
     * @param x - the argument of the polynomial.
     * 
     * @returns The value of the polynomial Lα,n for the argument x.
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

    /**
     * @brief Calculates the normalizing factor for spherical harmonics (Yl,m) or,
     * alternatively, for the associated Legendre polynomials (Pm,l) used in the definition of Yl,m.
     *
     * @tparam T - floating point type of the argument.
     * 
     * @param l - the degree of the polynomial (0 <= l),
     * @param m - the order of the polynomial (0 <= m <= l),
     *
     * @returns Normalizing factor according to the formula: `sqrt[ (2·l+1)/4π · (l-m)!/(l+m)! ]`.
     *
     * @exception std::invalid_argument - when the arguments (l,m) do not meet the prerequisites.
     * 
     */
    template<typename T>
    T ylm_norm_factor(unsigned l, unsigned m)
    {
        static_assert(
            std::is_floating_point<T>::value,
            "T ylm_norm_factor(unsigned l, unsigned m): template requires a floating-point type parameter T."
        );

        static constexpr T PI = mangh::math::PI_V<T>;

        if (m > l)
            throw std::invalid_argument(
                "T ylm_norm_factor(unsigned l, unsigned m): arguments (l, m) do not meet the condition: l >= m."
            );

        T prod = 1.0;
        for (unsigned j = l - m + 1; j <= l + m; ++j)
            prod *= j;

        return std::sqrt((2.0 * l + 1) / 4.0 / PI / prod);
    }

    /**
     * @brief Associated Legendre polynomial renormalized for use with spherical harmonics
     * (i.e., a polynomial with the normalization factor of `√[ (2·l+1)/4π · (l-m)!/(l+m)! ]`
     * and also the Condon-Shortley factor of `(-1)^m)`.
     *                                      
     * @tparam T - floating point type of the argument.
     * 
     * @param l - the degree of the polynomial (0 <= l),
     * @param m - the order of the polynomial (0 <= m <= l),
     * @param x - the argument (-1 <= x <= 1).
     * 
     * @returns The value of the renormalized Legendre polynomial
     * of degree `l` and order `m` for the argument `x`.
     * 
     * @exception std::invalid_argument - when arguments do not meet the prerequisites.
     * 
     * @remarks The following is a `plegendre` algorithm (from the book
     * Numerical Recipes in C++, chapter 6.7 Spherical Harmonics)
     * that implements equation (6.7.8):
     * 
     * `P~m,l = √[ (2·l+1)/4π · (l-m)!/(l+m)! ] · Pm,l`
     * 
     * where Pm,l are ordinary Legendre polynomials defined by equation (6.7.4):
     * 
     * `Pm,l(x) = (-1)^m · (1 - x²)^(m/2) · d^m/dx^m Pl(x)`.
     * 
     */
    template<typename T>
    T ylm_assoc_legendre(unsigned l, unsigned m, T x)
    {
        static_assert(
            std::is_floating_point<T>::value,
            "T ylm_assoc_legendre(unsigned l, unsigned m, T x): function template requires a floating-point type parameter T."
        );

        if (m > l)
            throw std::invalid_argument(
            "T ylm_assoc_legendre(unsigned l, unsigned m, T x): arguments (l, m) do not meet the condition: l >= m."
            );

        if ((x < -1) || (x > 1))
            throw std::invalid_argument(
            "T ylm_assoc_legendre(unsigned l, unsigned m, T x): argument x does not meet the condition: -1 <= x <= 1."
            );

        static constexpr T PI = mangh::math::PI_V<T>;

        // compute P~m,m
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

        // compute P~m,m+1
        T oldfact = std::sqrt(2.0 * m + 3.0);
        T Pmmp1 = x * oldfact * Pmm;
        if (l == m+1)
            return Pmmp1;

        // compute P~m,l (m+1 < l)
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
     * @brief Associated Legendre polynomial Pm,l(x) (orthonormal).
     * 
     * @tparam T - floating point type of the argument.
     * 
     * @param l - the degree of the polynomial (0 <= l),
     * @param m - the order of the polynomial (0 <= m <= l)
     * @param x - the argument (-1 <= x <= 1).
     * 
     * @exception std::invalid_argument - when arguments do not meet the prerequisites.
     * 
     * @returns Value of the associated Legendre polynomial (i.e. Pm,l(x)).
     * 
     * @remarks The function may overflow for m >= 80, or even sooner if l >> m.
     * @remarks See: Numerical Recipes in C++, 6.7. Spherical Harmonics, equation (6.7.4), algorithm "plgndr".
     * 
     */
    template<typename T>
    T assoc_legendre(unsigned l, unsigned m, T x)
    {
        static_assert(
            std::is_floating_point<T>::value,
            "T assoc_legendre(unsigned l, unsigned m, T x): template requires a floating-point type parameter T."
        );

        return ylm_assoc_legendre(l, m, x) / ylm_norm_factor<T>(l, m);
    }
}

#endif /* !MANGH_POLYNOMIALS_HPP */