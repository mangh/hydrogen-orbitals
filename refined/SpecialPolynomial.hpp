#ifndef SPECIAL_POLYNOMIALS_HPP
#define SPECIAL_POLYNOMIALS_HPP

#include <cassert>
#include <cmath>
#include <type_traits>
#include <stdexcept>

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
    T factorial(size_t n)
    {
        static_assert(std::is_floating_point<T>::value, "factorial(n) requires a floating-point type parameter.");
        return std::exp(std::lgamma(n + 1));
    }

    // constexpr inline size_t binom(size_t n, size_t k) noexcept
    // {
    //     return
    //         /* out of range? */ (k > n) ? 0 :           
    //         /* edge? */         (k == 0 || k == n) ? 1 :
    //         /* first? */        (k == 1 || k == n-1) ? n :
    //         /* recursive */     (k + k < n) ?
    //                                 (binom(n - 1, k - 1) * n) / k : // path to k = 1 is faster
    //                                 (binom(n - 1, k) * n) / (n - k);    // path to k = n-1 is faster
    // }

    // constexpr inline size_t binom2(size_t n, size_t k) noexcept
    // {
    //     return
    //         (k > n) ? 0 :
    //         (k == 0 || k == n) ? 1 :
    //         (k == 1 || k == n - 1) ? n :
    //             binom(n - 1, k - 1) * n / k;
    // }

    /**
     * @brief Binomial coefficient "n choose k" (n/k).
     * 
     * @param n - numerator argument,
     * @param k - denominator argument.
     * @returns The coefficient (n/k) in the expansion of the polynomial (x+y)^n.
     * 
     * @remarks See:
     * https://en.wikipedia.org/wiki/Binomial_coefficient
     * https://stackoverflow.com/questions/44718971/calculate-binomial-coffeficient-very-reliably
     * 
     */
    template<typename T>
    T binom(size_t n, size_t k)
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
                   * binom<double>(n, k)
                   * binom<double>(2 * n - 2 * k, n);
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
                * factorial<double>(n - (-m))
                / factorial<double>(n + (-m))
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

    /**
     * @brief Renormalized associated Legendre polynomial Pml(x)
     * 
     * @param l 
     * @param m 
     * @param x 
     * @returns double 
     * 
     * @pre 0 <= m <= l
     * @pre -1 <= x <= 1
     * 
     * @remark See: Numerical Recipes in C++. 6.7. Spherical Harmonics
     * 
     */
    template<typename T>
    T plegendre(int l, int m, T x)
    {
        static constexpr double PI = 3.141592653589793238462643383279502884L;

        if (m < 0 || m > l || abs(x) > 1.0)
            throw std::invalid_argument("plegendre: bad arguments.");

        double pmm = 1.0;
        if (m > 0)
        {
            double omx2 = (1.0 - x) * (1.0 + x);
            double fact = 1.0;
            for (int i = 1; i <= m; i++)
            {
                pmm *= omx2 * fact / (fact + 1.0);
                fact += 2.0;
            }
        }
        
        pmm = sqrt((2 * m + 1) * pmm / (4.0 * PI));

        if (m & 1)
            pmm = -pmm;

        if (l == m)
            return pmm;
        else
        {
            double pmmp1 = x * sqrt(2.0 * m + 3.0) * pmm;
            if (l == m+1)
                return pmmp1;
            else
            {
                double oldfact = sqrt(2.0 * m + 3.0);
                double pll = 0.0;   // ????
                for (int ll = m + 2; ll <= l; ll++)
                {
                    double fact = sqrt((4.0 * ll * ll - 1.0) / (ll * ll - m * m));
                    pll = (x * pmmp1 - pmm / oldfact) * fact;
                    oldfact = fact;
                    pmm = pmmp1;
                    pmmp1 = pll;
                }
                return pll;
            }
        }
    }

    /**
     * @brief Associated Legendre polynomial Pml(x) with standard normalization.
     * 
     * @param l 
     * @param m 
     * @param x 
     * @returns double 
     * 
     * @pre 0 <= m <= l
     * @pre -1 <= x <= 1
     * 
     * @remarks Overflow for m >= 80, or even sooner if l >> m.
     * @remarks See: Numerical Recipes in C++. 6.7. Spherical Harmonics
     * 
     */
    template<typename T>
    T plgndr(int l, int m, T x)
    {
        static constexpr double PI = 3.141592653589793238462643383279502884L;

        if (m < 0 || m > l || abs(x) > 1.0)
            throw std::invalid_argument("plgndr: bad arguments.");

        double prod = 1.0;
        for (int j = l - m + 1; j <= l + m; j++)
        {
            prod *= j;
        } 
        return sqrt(4.0 * PI * prod / (2 * l + 1)) * plegendre(l,m,x);
    }

}

#endif /* !SPECIAL_POLYNOMIALS_HPP */