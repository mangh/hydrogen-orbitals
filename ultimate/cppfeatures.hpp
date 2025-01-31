#ifndef CPP_FEATURES_H
#define CPP_FEATURES_H

/* 
 * __cplusplus
 *
 *  Macro constant: denotes the version of C++ standard that is being used, expands to value:
 * 
 *  199711L(until C++11),
 *  201103L(C++11),
 *  201402L(C++14),
 *  201703L(C++17),
 *  202002L(C++20), or
 *  202302L(C++23)
 * 
 *  See: https://en.cppreference.com/w/cpp/preprocessor/replace
 *
 *
 * _MSVC_LANG
 *
 *  Defined as an integer literal that specifies the C++ language standard targeted by the compiler
 *
 *  See: https://learn.microsoft.com/en-us/cpp/preprocessor/predefined-macros?view=msvc-170
 * 
 */

#ifdef __cplusplus

#  ifndef CPP_VERSION
#    if (__cplusplus >= 202302L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 202302L)
#      define CPP_VERSION 2023
#    elif (__cplusplus >= 202002L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)
#      define CPP_VERSION 2020
#    elif (__cplusplus >= 201703L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L)
#      define CPP_VERSION 2017
#    elif (__cplusplus >= 201402L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 201402L)
#      define CPP_VERSION 2014
#    elif (__cplusplus >= 201103L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 201103L)
#      define CPP_VERSION 2011
#    else
#      define CPP_VERSION 1997
#    endif
#  endif // CPP_VERSION

#  if CPP_VERSION > 2011
#    define CONSTEXPR_SINCE_CXX14 constexpr
#  else
#    define CONSTEXPR_SINCE_CXX14
#  endif

#  if CPP_VERSION > 2014
#    define CONSTEXPR_SINCE_CXX17 constexpr
#  else
#    define CONSTEXPR_SINCE_CXX17
#  endif

#  if CPP_VERSION > 2017
#    define CONSTEXPR_SINCE_CXX20 constexpr
#  else
#    define CONSTEXPR_SINCE_CXX20
#  endif

#  if CPP_VERSION > 2020
#    define CONSTEXPR_SINCE_CXX23 constexpr
#  else
#    define CONSTEXPR_SINCE_CXX23
#  endif

#endif // __cplusplus

#endif  /* !CPP_FEATURES_H */