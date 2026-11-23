#ifndef MANGH_MATH_ASSOC_LEGENDRE_HPP
#define MANGH_MATH_ASSOC_LEGENDRE_HPP

#include <stdexcept>
#include <type_traits>
#include <cmath>

#include "constants.hpp"

namespace mangh { namespace math
{
    /**
     * @brief Calculates the normalizing factor for spherical harmonics `Yl,m` or,
     * alternatively, for the associated Legendre polynomials `Pl,m` used in the definition of `Yl,m`.
     *
     * @tparam T - floating point type of the return value.
     * 
     * @param l - the degree of the polynomial (0 <= l),
     * @param m - the order of the polynomial (0 <= m <= l),
     *
     * @returns `Yl,m` normalizing factor = `√[ (2·l+1)/4π · (l-m)!/(l+m)! ]`.
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
     * @brief Associated Legendre polynomial `P~l,m` renormalized for use with spherical harmonics (i.e.,
     * a polynomial with the normalizing factor of `√[ (2·l+1)/4π · (l-m)!/(l+m)! ]` and the Condon-Shortley factor of `(-1)^m)`.
     *                                      
     * @tparam T - floating point type of the argument.
     * 
     * @param l - the degree of the polynomial (0 <= l),
     * @param m - the order of the polynomial (0 <= m <= l),
     * @param x - the argument (-1 <= x <= 1).
     * 
     * @returns The value of the renormalized Legendre `P~l,m` at `x`.
     * 
     * @exception std::invalid_argument - when arguments do not meet the prerequisites.
     * 
     * @remarks The function is a copy of the `plegendre` algorithm
     * (from the book Numerical Recipes in C++, chapter 6.7 Spherical Harmonics)
     * that implements equation (6.7.8):
     * 
     * `P~l,m = √[ (2·l+1)/4π · (l-m)!/(l+m)! ] · Pl,m` 
     * 
     * where `Pl,m` is an ordinary Legendre polynomials defined by equation (6.7.4):
     * 
     * `Pl,m(x) = (-1)^m · (1 - x²)^(m/2) · d^m/dx^m Pl(x)`.
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
 
         if (m % 2)
             Pmm = -Pmm;
 
         if (l == m)
             return Pmm;
 
         // compute P~m,m+1
         T oldfact = std::sqrt(2.0 * m + 3.0);
         T Pmmp1 = x * oldfact * Pmm;
         if (l == m+1)
             return Pmmp1;
 
         // compute P~l,m (m+1 < l)
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
      * @brief Associated Legendre polynomial `Pl,m(x)` (unnormalized, with the Condon-Shortley phase `(-1)^m`).
      * 
      * @tparam T - floating point type of argument and return value.
      * 
      * @param l - the degree of the polynomial (0 <= l),
      * @param m - the order of the polynomial (0 <= m <= l)
      * @param x - the argument (-1 <= x <= 1).
      * 
      * @exception std::invalid_argument - when arguments do not meet the prerequisites.
      * 
      * @returns Value of the associated Legendre polynomial `Pl,m` at `x`.
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
}}
    
#endif /* !MANGH_MATH_ASSOC_LEGENDRE_HPP */
