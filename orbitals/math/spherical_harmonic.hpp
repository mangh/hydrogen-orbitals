#ifndef MANGH_MATH_SPHERICAL_HARMONIC_HPP
#define MANGH_MATH_SPHERICAL_HARMONIC_HPP

#include <cmath>
#include <complex>

#include "assoc_legendre.hpp"

namespace mangh { namespace math
{
    /**
     * @brief Spherical harmonic `Yl,m(θ,ϕ) = P~l,m(cosθ) · e^(imϕ)`.
     *
     * @param l - the degree of the polynomial (0 <= l),
     * @param m - the order of the polynomial (0 <= m <= l),
     * @param theta - polar angle 0 ≤ θ ≤ π (from the positive Z axis),
     * @param phi - azimuthal angle 0 ≤ ϕ < 2π (from the X-axis, in the XY-plane).
     *
     * @returns The value of the spherical harmonic function `Yl,m` for the arguments (θ,ϕ).
     * 
     * @remarks The algorithm uses associated Legendre polynomials `P~l,m`,
     * @remarks which include both the Condon-Shortley factor and the normalizing
     * @remarks factor for spherical harmonics, i.e. the product:
     * @remarks 
     * @remarks     `(-1)^m · sqrt[ (2·l+1)/4π · (l-m)!/(l+m)! ]`.
     * @remarks 
     * @remarks It uses also the relationship:
     * @remarks 
     * @remarks     `Yl,-m(θ,ϕ) = (-1)^m · conj(Yl,m(θ,ϕ))`
     * @remarks 
     * @remarks See: Numerical Recipes in C++, Chapter 6.7, Spherical Harmonics.
     * 
     */
    template<typename T>
    std::complex<T> spherical_harmonic(unsigned l, int m, T theta, T phi)
    {
        unsigned M = m < 0 ? -m : m;
        std::complex<T> phase{ std::cos(M * phi), std::sin(M * phi) };
        return
            ylm_assoc_legendre(l, M, std::cos(theta)) *
                (m >= 0 ? phase :
                    (m % 2 ? -std::conj(phase) : 
                        std::conj(phase)));
    }

    /**
     * @brief Real-valued spherical harmonic basis function `Yl,m(θ,ϕ)`.
     *
     * @param l - the degree of the polynomial (0 <= l),
     * @param m - the order of the polynomial (-l <= m <= l),
     * @param theta - polar angle 0 ≤ θ ≤ π (from the positive Z axis),
     * @param phi - azimuthal angle 0 ≤ ϕ < 2π (from the X-axis, in the XY-plane).
     *
     * @returns The value of the real-valued spherical harmonic basis function `Yl,m` for the arguments (θ,ϕ).
     * 
     */
    template <typename T>
    T spherical_real_basis(unsigned l, int m, T theta, T phi)
    {
        unsigned M = m < 0 ? -m : m;
        auto ylm = ylm_assoc_legendre(l, M, std::cos(theta));
        return 
            (m < 0) ? std::sqrt(2.) * (m % 2 ? -ylm : ylm) * std::sin(M * phi) :
                (m > 0) ? std::sqrt(2.) * (m % 2 ? -ylm : ylm) * std::cos(M * phi) :
                    ylm;
    }
}}

#endif /* !MANGH_MATH_SPHERICAL_HARMONIC_HPP */
