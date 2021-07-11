/**
 *    ______             __     __  _____ _____
 *   / __/ /____ _____  / /__  /  |/  / // / _ \
 *  _\ \/  '_/ // / _ \/  '_/ / /|_/ / _  / // /
 * /___/_/\_\\_,_/_//_/_/\_\ /_/  /_/_//_/____/
 *
 * Copyright (c) 2019 Oleg Butakov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once
#ifndef FEATHERS_BASE_HH_
#define FEATHERS_BASE_HH_

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES 1
#endif
#ifndef NOMINMAX
#define NOMINMAX 1
#endif

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cassert>

#include <array>
#include <vector>
#include <string>
#include <memory>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <type_traits>

//#define GLM_FORCE_CTOR_INIT 1
#include <glm/glm.hpp>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** C++03 support. */
/** @{ */
#if (__cplusplus >= 199711L) || (_MSVC_LANG >= 199711L)
#define FEATHERS_HAS_CPP03 1
#else
#define FEATHERS_HAS_CPP03 0
#endif
/** @} */

/** C++11 support. */
/** @{ */
#if (__cplusplus >= 201103L) || (_MSVC_LANG >= 201103L)
#define FEATHERS_HAS_CPP11 1
#else
#define FEATHERS_HAS_CPP11 0
#endif
/** @} */

/** C++14 support. */
/** @{ */
#if (__cplusplus >= 201402L) || (_MSVC_LANG >= 201402L)
#define FEATHERS_HAS_CPP14 1
#else
#define FEATHERS_HAS_CPP14 0
#endif
/** @} */

/** C++17 support. */
/** @{ */
#if (__cplusplus >= 201703L) || (_MSVC_LANG >= 201703L)
#define FEATHERS_HAS_CPP17 1
#else
#define FEATHERS_HAS_CPP17 0
#endif
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** OpenMP 2.0 support. */
/** @{ */
#if (defined(_OPENMP) && (_OPENMP >= 200203))
#define FEATHERS_HAS_OPENMP2 1
#else
#define FEATHERS_HAS_OPENMP2 1
#endif
/** @} */

/** OpenMP 3.0 support. */
/** @{ */
#if (defined(_OPENMP) && (_OPENMP >= 200805))
#define FEATHERS_HAS_OPENMP3 1
#else
#define FEATHERS_HAS_OPENMP3 1
#endif
/** @} */

/** OpenMP 4.0 support. */
/** @{ */
#if (defined(_OPENMP) && (_OPENMP >= 201307))
#define FEATHERS_HAS_OPENMP4 1
#else
#define FEATHERS_HAS_OPENMP4 1
#endif
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** Convert token to string. */
/** @{ */
#define FEATHERS_TO_STRING_(x) #x
#define FEATHERS_TO_STRING(x) FEATHERS_TO_STRING_(x)
/** @} */

/** Concatenate tokens. */
/** @{ */
#define FEATHERS_CONCAT_(x, y) x##y
#define FEATHERS_CONCAT(x, y) FEATHERS_CONCAT_(x, y)
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#define FEATHERS_NOT_USED(x)

/** Pragma. */
/** @{ */
#ifdef _MSC_VER
#define FEATHERS_PRAGMA(...) __pragma(__VA_ARGS__)
#else
#define FEATHERS_PRAGMA(...) _Pragma(FEATHERS_TO_STRING(__VA_ARGS__))
#endif
/** @} */

/** A bit more aggressive inline macro. */
/** @{ */
#ifdef _MSC_VER
#define SKUNK_INLINE inline __forceinline
#else
#define SKUNK_INLINE inline __attribute__((always_inline))
#endif
/** @} */

/** Assume macro: an optimization hint for the compiler. */
/** @{ */
#ifdef _MSC_VER
#define FEATHERS_ASSUME(x) __assume(x)
#else
#define FEATHERS_ASSUME(x) do { if (!(x)) { __builtin_unreachable(); } } while (false)
#endif
/** @} */

/** Likely/Unlikely macro: a branching optimization hint for the compiler. */
/** @{ */
#ifdef _MSC_VER
#define FEATHERS_LIKELY(x) (x)
#define FEATHERS_UNLIKELY(x) (x)
#else
#define FEATHERS_LIKELY(x) __builtin_expect((x), 1)
#define FEATHERS_UNLIKELY(x) __builtin_expect((x), 0)
#endif
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** Compile-time expression macro. */
/** @{ */
#if FEATHERS_HAS_CPP14
#define SKUNK_CONSTEXPR constexpr
#else
#define SKUNK_CONSTEXPR const
#endif
/** @} */

/** Compile-time if statement macro. */
/** @{ */
#if FEATHERS_HAS_CPP17
#define SKUNK_IF_CONSTEXPR if constexpr
#else
#define SKUNK_IF_CONSTEXPR if
#endif
/** @} */

/** Compile-time assertion macro.
 ** Assertion message is optional. */
/** @{ */
#if FEATHERS_HAS_CPP17
#define SKUNK_STATIC_ASSERT(x) static_assert(x)
#else
#define SKUNK_STATIC_ASSERT(x) static_assert(x, FEATHERS_TO_STRING(x))
#endif
/** @} */

/**************************************************************************/
/**************************************************************************/

#define FEATHERS_ENSURE(x) do { if(!(x)) { \
        std::fprintf(stderr, "\nAssertion failed:\n%s:%d %s: \"%s\".\n", \
                     __FILE__, __LINE__, __FUNCTION__, FEATHERS_TO_STRING(x)); \
        std::abort(); \
    } } while (false)

/** Assertion macro. */
/** @{ */
#ifdef NDEBUG
#define FEATHERS_ASSERT(x) FEATHERS_ASSUME(x)
#else
#define FEATHERS_ASSERT(x) FEATHERS_ENSURE(x)
#endif
/** @} */

/** Fatal assertion macro. */
/** @{ */
#ifdef NDEBUG
#define FEATHERS_ASSERT_FALSE(m) do { std::abort(); } while (false)
#else
#define FEATHERS_ASSERT_FALSE(m) do { \
        std::fprintf(stderr, "\nFatal assertion failed:\n%s:%d %s: %s", \
                     __FILE__, __LINE__, __FUNCTION__, m); \
        std::abort(); \
    } while (false)
#endif
/** @} */

/** Not implemented macro. */
#define FEATHERS_NOT_IMPLEMENTED() FEATHERS_ASSERT_FALSE("not implemented.")
/** Not implemented macro. */
#define FEATHERS_NOT_REACHABLE() FEATHERS_ASSERT_FALSE("not reachable.")

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Unsigned byte type. */
using byte_t = std::uint8_t;

/** A regular size type. */
using size_t = std::size_t;
/** A regular pointer difference type. */
using ptrdiff_t = std::ptrdiff_t;

/** Signed integer type, 32-bit wide. */
using int_t = std::int32_t;
/** Unsigned integer type, 32-bit wide. */
using uint_t = std::uint32_t;

/** Invalid index. */
static constexpr uint_t npos = std::numeric_limits<uint_t>::max();

/** Check if index is npos. */
static constexpr auto is_npos(uint_t ind) {
    return ind == npos;
}   // is_npos
/** Check if index is not npos. */
static constexpr auto is_not_npos(uint_t ind) {
    return ind != npos;
}   // is_not_npos

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#ifndef SKUNK_CONFIG_REAL_IS_DOUBLE
#define SKUNK_CONFIG_REAL_IS_DOUBLE 1
#endif

/** Floating-point type, vector and matrix types. */
/** @{ */
#if SKUNK_CONFIG_REAL_IS_DOUBLE
using real_t = double;
using vec2_t = glm::dvec2; using vec3_t = glm::dvec3; using vec4_t = glm::dvec4;
using mat2_t = glm::dmat2; using mat3_t = glm::dmat3; using mat4_t = glm::dmat4;
#else
using real_t = float;
using vec2_t = glm::vec2; using vec3_t = glm::vec3; using vec4_t = glm::vec4;
using mat2_t = glm::mat2; using mat3_t = glm::mat3; using mat4_t = glm::mat4;
#endif
/** @} */

/** A @f$e@f$ constant. */
/** @{ */
#ifdef M_E
static constexpr real_t c_e(M_E);
#else
static constexpr real_t c_e(std::exp(1.0));
#endif
/** @} */
/** A @f$\log_2(e)@f$ constant. */
/** @{ */
#ifdef M_LOG2E
static constexpr real_t c_log2e(M_LOG2E);
#else
static constexpr real_t c_log2e(std::log2(m_e));
#endif
/** @} */
/** A @f$\log_10(e)@f$ constant. */
/** @{ */
#ifdef M_LOG10E
static constexpr real_t c_log10e(M_LOG10E);
#else
static constexpr real_t c_log10e(std::log10(m_e));
#endif
/** @} */
/** A @f$\log(2)@f$ constant. */
/** @{ */
#ifdef M_LN2
static constexpr real_t c_ln2(M_LN2);
#else
static constexpr real_t c_ln2(std::log(2.0));
#endif
/** @} */
/** A @f$\log(10)@f$ constant. */
/** @{ */
#ifdef M_LN10
static constexpr real_t c_ln10(M_LN10);
#else
static constexpr real_t c_ln10(std::log(10.0));
#endif
/** @} */

/** A @f$\pi@f$ constant. */
/** @{ */
#ifdef M_PI
static constexpr real_t c_pi(M_PI);
#else
static constexpr real_t c_pi(4.0*std::atan(1.0));
#endif
/** @} */
/** A @f$\frac{\pi}{2}@f$ constant. */
/** @{ */
#ifdef M_PI_2
static constexpr real_t c_pi_2(M_PI_2);
#else
static constexpr real_t c_pi_2(2.0*std::atan(1.0));
#endif
/** A @f$\frac{\pi}{4}@f$ constant. */
/** @{ */
#ifdef M_PI_4
static constexpr real_t c_pi_4(M_PI_4);
#else
static constexpr real_t c_pi_4(1.0*std::atan(1.0));
#endif
/** @} */
/** A @f$\frac{1}{\pi}@f$ constant. */
/** @{ */
#ifdef M_1_PI
static constexpr real_t c_1_pi(M_1_PI);
#else
static constexpr real_t c_1_pi(1.0/m_pi);
#endif
/** @} */
/** A @f$\frac{2}{\pi}@f$ constant. */
/** @{ */
#ifdef M_2_PI
static constexpr real_t c_2_pi(M_2_PI);
#else
static constexpr real_t c_2_pi(1.0/m_pi_2);
#endif
/** @} */
/** A @f$\frac{2}{\sqrt{\pi}}@f$ constant. */
/** @{ */
#ifdef M_2_SQRTPI
static constexpr real_t c_2_sqrtpi(M_2_SQRTPI);
#else
static constexpr real_t c_2_sqrtpi(2.0/std::sqrt(m_pi));
#endif
/** @} */

/** A @f$\sqrt{2}@f$ constant. */
/** @{ */
#ifdef M_SQRT2
static constexpr real_t c_sqrt2(M_SQRT2);
#else
static constexpr real_t c_sqrt2(std::sqrt(2.0));
#endif
/** @} */
/** A @f$\frac{1}{\sqrt{2}}@f$ constant. */
/** @{ */
#ifdef M_SQRT1_2
static constexpr real_t c_sqrt1_2(M_SQRT1_2);
#else
static constexpr real_t c_sqrt1_2(std::sqrt(0.5));
#endif
/** @} */

/** Compute pseudo-inverse value. */
template<typename type_t>
static constexpr type_t safe_inverse(type_t x) {
    constexpr type_t zero{};
    return x == zero ? zero : type_t(1)/x;
}   // safe_inverse

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Simple shortcut for @c std::enable_shared_from_this. */
template<typename type_t>
using tObject = std::enable_shared_from_this<type_t>;

template<typename obj_t, typename obj_ptr_t>
SKUNK_INLINE auto shared_from_this(const obj_ptr_t& obj) {
    return std::static_pointer_cast<obj_t>(obj->shared_from_this());
}   // shared_from_this

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/* TODO: Remove me. */
using int_t = feathers::int_t;
using uint_t = feathers::uint_t;
using real_t = feathers::real_t;

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
