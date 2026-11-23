#ifndef MANGH_MATH_KAHAN_SUM_HPP
#define MANGH_MATH_KAHAN_SUM_HPP

#include <type_traits>

namespace mangh { namespace math
{
    /**
     * @brief Compensated (Kahan) summation: keeps the running rounding error
     * and feeds it back into the next addition, so that the total does not drift
     * when many small terms are added to a large sum.
     *
     * @tparam T - floating point type of the terms.
     *
     * @remarks Usage: `KahanSum<double> s; for (...) s += term; double total = s.value();`
     */
    template<typename T>
    class KahanSum
    {
        static_assert(std::is_floating_point_v<T>, "KahanSum<T>: T must be a floating-point type.");

    public:

        constexpr KahanSum& operator+=(T addend) noexcept
        {
            const T y = addend - m_error;
            const T t = m_total + y;
            m_error = (t - m_total) - y;
            m_total = t;
            return *this;
        }

        /// @brief The sum of the terms added so far.
        constexpr T value() const noexcept { return m_total; }

    private:
        T m_total{};
        T m_error{};
    };
}}

#endif /* !MANGH_MATH_KAHAN_SUM_HPP */
