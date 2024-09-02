#ifndef POLYNOMIALS_HPP
#define POLYNOMIALS_HPP

#include <cassert>
#include <cmath>
#include <type_traits>

namespace SpecialPolynomial
{
    template<typename N>
    N factorial(N n)
    {
        static_assert(std::is_integral<N>::value, "factorial() function requires an integer type parameter.");
        return (n <= 0) ? 1 : n * factorial(n - 1);
    }

    /**
     * @brief Legendre polynomial Pl(x).
     * 
     * @tparam T 
     * @param x 
     * @param l 
     * @returns T 
     * @remarks Recursive definition:
     *      
     *  1    l/2+1  (-1)^k∙(2l-2k)!∙x^(l-2k)
     * --- ∙   Σ    ------------------------
     * 2^l    k=0      k!∙(l-k)!∙(l-2k)!
     *      
     * See: https://mathworld.wolfram.com/LegendrePolynomial.html.
     *      
     */
    template<typename T>
    T Pl(T x, int l)
    {
        static_assert(std::is_floating_point<T>::value, "Legendre polynomial Pl(x,l) requires a floating-point type parameter x.");
        T total = 0.0;
        for (int k = 0; k < l / 2; ++k)
        {
            total += pow(-1.0, k)
                   * pow(x, l - 2 * k)
                   * factorial(2 * l - 2 * k)
                   / factorial(k)
                   / factorial(l - k)
                   / factorial(l - 2 * k);
        }
        return total / pow(2.0, l);
    }

    template<typename T>
    T dm_Pl_dxm(T x, int m, int l)
    {
        assert(0 <= m);
        T dx = 0.0001;
        return
            (m == 0) ? Pl(x, l) :
            (m == 1) ? (Pl(x + dx, l) - Pl(x, l)) / dx :
            (dm_Pl_dxm(x + dx, m - 1, l) - dm_Pl_dxm(x, m - 1, l)) / dx;
    }

    /**
     * @brief Associated Legendre polynomial Pml(m,l,x).
     * 
     * @tparam T 
     * @param m 
     * @param l 
     * @param x 
     * @returns T 
     */
    template<typename T>
    T Pml(int m, int l, T x)
    {
        static_assert(std::is_floating_point<T>::value, "Associated Legendre polynomial Pml(m,l,x) requires a floating-point type parameter x.");
        return(m >= 0) ?
            pow(-1.0, m) * pow(1.0 - x * x, m / 2.0) * dm_Pl_dxm(x, m, l) :
            pow(-1.0, (-m)) * factorial(l - (-m)) / factorial(l + (-m)) * Pml((-m), l, x);
    }

    /**
     * @brief Generalized (or associated) Laguerre polynomial.
     * 
     * @tparam T 
     * @param x 
     * @param a 
     * @param k 
     * @returns T 
     * 
     * @remarks Recursive definition:
     * (n+1) ∙ L[n+1](x) = (2n+1-x) ∙ L[n](x) - n ∙ L[n-1](x).
     */
    template<typename T>
    T Ln(T x, int a, int k)
    {
        static_assert(std::is_floating_point<T>::value, "Laguerre polynomial requires a floating-point type parameter x.");
        return
            (k == 0) ? 1.0 :
            (k == 1) ? 1.0 + a - x :
            ((2.0 * (k - 1) + 1.0 + a - x) * Ln(x, a, k - 1) - (k - 1 + a) * Ln(x, a, k - 2)) / k;
    }
}

#endif /* !POLYNOMIALS_HPP */