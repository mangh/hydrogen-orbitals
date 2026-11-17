#ifndef MANGH_MATH_CONSTANTS_HPP
#define MANGH_MATH_CONSTANTS_HPP

#include <type_traits>

namespace mangh { namespace math
{
    ///////////////////////////////////////////////////////////////////
    ///
    ///     Constants
    ///

    /**
     * @brief Expression to check use of floating-point argument
     * (after /usr/include/c++/12/numbers in Debian GCC).
     * 
     * @tparam _Tp - the type of the argument to check (whether it is a floating point type).
     */
    template<typename _Tp>
    using _Enable_if_floating = std::enable_if_t<std::is_floating_point_v<_Tp>, _Tp>;

    /**
     * 
     * @brief PI constant (π)
     *
     * @remark OS         | compiler | header file                        | PI constant
     * @remark ---------- + -------- + ---------------------------------- + ---------------------------------------
     * @remark Windows 11 | msvc     | %VCToolsInstallDir%include\numbers | 3.141592653589793
     * @remark Windows 11 | msvc     | %VCToolsInstallDir%include\random  | 3.14159265358979323846264338327950288L
     * @remark Windows 11 | msvc     | %VCToolsInstallDir%include\complex | 3.1415926535897932384626433832795029L
     * @remark Debian     | gcc      | /usr/include/c++/12/numbers        | 3.141592653589793238462643383279502884L
     * @remark FreeBSD    | clang    | /usr/include/c++/v1/numbers        | 3.141592653589793238462643383279502
     */ 
    template<typename _Tp>
    inline constexpr _Tp PI_V
        = _Enable_if_floating<_Tp>(3.141592653589793238462643383279502884L);

}}

#endif /* !MANGH_MATH_CONSTANTS_HPP */
