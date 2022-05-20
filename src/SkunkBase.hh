/*
 *  ______  ______   ______   ______  __  __   ______   ______   ______
 * /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\
 * \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \
 *  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\
 *   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
 *
 * Copyright (c) 2021 Oleg Butakov
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

#include <glm/glm.hpp>

#if 0
#include <ranges>
namespace ranges = std::ranges;
namespace views = std::views;
#else
#include <range/v3/all.hpp>
namespace views = ranges::views;
#endif

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** C++03 support. */
/** @{ */
#if (__cplusplus >= 199711L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 199711L)
#define FEATHERS_HAS_CPP_03 1
#else
#define FEATHERS_HAS_CPP_03 0
#endif
/** @} */

/** C++11 support. */
/** @{ */
#if (__cplusplus >= 201103L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 201103L)
#define FEATHERS_HAS_CPP_11 1
#else
#define FEATHERS_HAS_CPP_11 0
#endif
/** @} */

/** C++14 support. */
/** @{ */
#if (__cplusplus >= 201402L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 201402L)
#define FEATHERS_HAS_CPP_14 1
#else
#define FEATHERS_HAS_CPP_14 0
#endif
/** @} */

/** C++17 support. */
/** @{ */
#if (__cplusplus >= 201703L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L)
#define FEATHERS_HAS_CPP_17 1
#else
#define FEATHERS_HAS_CPP_17 0
#endif
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/* Force OpenMP 2.0 */
#ifndef FEATHERS_CONFIG_FORCE_DISABLE_OPENMP
#define FEATHERS_CONFIG_FORCE_DISABLE_OPENMP 0
#endif

/* Force OpenMP 2.0 */
#ifndef FEATHERS_CONFIG_FORCE_OPENMP_2_0
#define FEATHERS_CONFIG_FORCE_OPENMP_2_0 0
#endif

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#ifndef FEATHERS_HAS_TBB
#define FEATHERS_HAS_TBB 0
#endif

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

// TODO:
#define FEATHERS_DOXYGEN 0

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
#if FEATHERS_HAS_CPP_14
#define SKUNK_CONSTEXPR constexpr
#else
#define SKUNK_CONSTEXPR const
#endif
/** @} */

/** Compile-time if statement macro. */
/** @{ */
#if FEATHERS_HAS_CPP_17
#define SKUNK_IF_CONSTEXPR if constexpr
#else
#define SKUNK_IF_CONSTEXPR if
#endif
/** @} */

/** Compile-time assertion macro.
 ** Assertion message is optional. */
/** @{ */
#if FEATHERS_HAS_CPP_17
#define SKUNK_STATIC_ASSERT(x) static_assert(x)
#else
#define SKUNK_STATIC_ASSERT(x) static_assert(x, FEATHERS_TO_STRING(x))
#endif
/** @} */

#define FEATHERS_DEPRECATED //[[deprecated("")]]

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#define FEATHERS_PASS(...) __VA_ARGS__

/** @{ */
#define FEATHERS_CONST_OVERLOAD_R_T(T, type, const_type, method_name, arguments, ...) \
    /** @{ */ \
    T type method_name arguments __VA_ARGS__ \
    T const_type method_name arguments const __VA_ARGS__ \
    /** @} */
#define FEATHERS_CONST_OVERLOAD_R(type, const_type, method_name, arguments, ...) \
    FEATHERS_CONST_OVERLOAD_R_T( \
        /*empty*/, type, const_type, method_name, arguments, __VA_ARGS__)
#define FEATHERS_CONST_OVERLOAD_T(T, type, method_name, arguments, ...) \
    FEATHERS_CONST_OVERLOAD_R_T( \
        T, type, const type, method_name, arguments, __VA_ARGS__)
#define ConstOverload(type, method_name, arguments, ...) \
    FEATHERS_CONST_OVERLOAD_T( \
        /*empty*/, type, method_name, arguments, __VA_ARGS__)
/** @} */

#define StormAutoConstOverload(name, arguments, ...) \
  /** @{ */ \
  auto name arguments __VA_ARGS__ \
  auto name arguments const __VA_ARGS__ \
  /** @} */

/* https://stackoverflow.com/a/67374211 */
#if (!defined(__PRETTY_FUNCTION__) && !defined(__GNUC__))
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#define FEATHERS_ENSURE(x) do { if(!(x)) { \
        std::fprintf( \
            stderr, "\nAssertion failed:\n%s:%d %s: \"%s\".\n", \
            __FILE__, __LINE__, __PRETTY_FUNCTION__, FEATHERS_TO_STRING(x)); \
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

#define StormAssert FEATHERS_ASSERT
#define StormEnsure FEATHERS_ENSURE

/** Fatal assertion macro. */
#define FEATHERS_ERROR_STOP(message) do { \
        std::fprintf( \
            stderr, "\nFatal assertion failed:\n%s:%d %s: %s", \
            __FILE__, __LINE__, __PRETTY_FUNCTION__, message); \
        std::abort(); \
    } while (false)

/** Not implemented macro. */
#define FEATHERS_NOT_IMPLEMENTED(...) \
    FEATHERS_ERROR_STOP("not implemented. " __VA_ARGS__)
/** Not implemented macro. */
#define FEATHERS_NOT_REACHABLE(...) \
    FEATHERS_ERROR_STOP("not reachable. " __VA_ARGS__)

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

namespace feathers {

/** Unsigned byte type. */
using byte_t = std::uint8_t;

/** A regular size type. */
using size_t = std::size_t;
/** A regular pointer difference type. */
using ptrdiff_t = std::ptrdiff_t;

/** Invalid index. */
static constexpr size_t npos = std::numeric_limits<size_t>::max();

/** Check if index is npos. */
template<class Value>
static constexpr bool is_npos(Value ind) {
  return ind == npos;
} // is_npos
/** Check if index is not npos. */
template<class Value>
static constexpr bool is_not_npos(Value ind) {
  return ind != npos;
} // is_not_npos

/** Min value functor. */
template<typename tValue>
class tMinFunc {
public:
  constexpr tValue operator()(tValue value_1, tValue value_2) const {
    return std::min(value_1, value_2);
  }
}; // class tMinFunc

/** Max value functor. */
template<typename tValue>
class tMaxFunc {
public:
  constexpr tValue operator()(tValue value_1, tValue value_2) const {
    return std::max(value_1, value_2);
  }
}; // class tMaxFunc

/** Min-max value functor. */
template<typename tValue>
class tMinMaxFunc {
public:
  constexpr auto operator()(tValue value_1, tValue value_2) const {
    return std::minmax(value_1, value_2);
  }
  constexpr auto operator()(tValue value_1,
                            const std::pair<tValue, tValue>& value_2) const {
    return std::make_pair(std::min(value_1, value_2.first),
                          std::max(value_1, value_2.second));
  }
  constexpr auto operator()(const std::pair<tValue, tValue>& value_1,
                            tValue value_2) const {
    return std::make_pair(std::min(value_1.first, value_2),
                          std::max(value_1.second, value_2));
  }
  constexpr auto operator()(const std::pair<tValue, tValue>& value_1,
                            const std::pair<tValue, tValue>& value_2) const {
    return std::make_pair(std::min(value_1.first, value_2.first),
                          std::max(value_1.second, value_2.second));
  }
}; // class tMinMaxFunc

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

/** Maximum real number constant. */
static constexpr real_t huge = std::numeric_limits<real_t>::max();

/** Quiet NaN constant. */
static constexpr real_t qnan = std::numeric_limits<real_t>::quiet_NaN();

/** A @f$e@f$ constant. */
/** @{ */
#ifdef M_E
static constexpr real_t c_e(M_E);
#else
static const real_t c_e(std::exp(1.0));
#endif
/** @} */
/** A @f$\log_2(e)@f$ constant. */
/** @{ */
#ifdef M_LOG2E
static constexpr real_t c_log2e(M_LOG2E);
#else
static const real_t c_log2e(std::log2(m_e));
#endif
/** @} */
/** A @f$\log_10(e)@f$ constant. */
/** @{ */
#ifdef M_LOG10E
static constexpr real_t c_log10e(M_LOG10E);
#else
static const real_t c_log10e(std::log10(m_e));
#endif
/** @} */
/** A @f$\log(2)@f$ constant. */
/** @{ */
#ifdef M_LN2
static constexpr real_t c_ln2(M_LN2);
#else
static const real_t c_ln2(std::log(2.0));
#endif
/** @} */
/** A @f$\log(10)@f$ constant. */
/** @{ */
#ifdef M_LN10
static constexpr real_t c_ln10(M_LN10);
#else
static const real_t c_ln10(std::log(10.0));
#endif
/** @} */

/** A @f$\pi@f$ constant. */
/** @{ */
#ifdef M_PI
static constexpr real_t c_pi(M_PI);
#else
static const real_t c_pi(4.0*std::atan(1.0));
#endif
/** @} */
/** A @f$\frac{\pi}{2}@f$ constant. */
/** @{ */
#ifdef M_PI_2
static constexpr real_t c_pi_2(M_PI_2);
#else
static const real_t c_pi_2(2.0*std::atan(1.0));
#endif
/** A @f$\frac{\pi}{4}@f$ constant. */
/** @{ */
#ifdef M_PI_4
static constexpr real_t c_pi_4(M_PI_4);
#else
static const real_t c_pi_4(1.0*std::atan(1.0));
#endif
/** @} */
/** A @f$\frac{1}{\pi}@f$ constant. */
/** @{ */
#ifdef M_1_PI
static constexpr real_t c_1_pi(M_1_PI);
#else
static const real_t c_1_pi(1.0/m_pi);
#endif
/** @} */
/** A @f$\frac{2}{\pi}@f$ constant. */
/** @{ */
#ifdef M_2_PI
static constexpr real_t c_2_pi(M_2_PI);
#else
static const real_t c_2_pi(1.0/m_pi_2);
#endif
/** @} */
/** A @f$\frac{2}{\sqrt{\pi}}@f$ constant. */
/** @{ */
#ifdef M_2_SQRTPI
static constexpr real_t c_2_sqrtpi(M_2_SQRTPI);
#else
static const real_t c_2_sqrtpi(2.0/std::sqrt(m_pi));
#endif
/** @} */

/** A @f$\sqrt{2}@f$ constant. */
/** @{ */
#ifdef M_SQRT2
static constexpr real_t c_sqrt2(M_SQRT2);
#else
static const real_t c_sqrt2(std::sqrt(2.0));
#endif
/** @} */
/** A @f$\frac{1}{\sqrt{2}}@f$ constant. */
/** @{ */
#ifdef M_SQRT1_2
static constexpr real_t c_sqrt1_2(M_SQRT1_2);
#else
static const real_t c_sqrt1_2(std::sqrt(0.5));
#endif
/** @} */

/** Compute pseudo-inverse value. */
template<typename tValue>
constexpr tValue safe_inverse(tValue x) {
  return x == tValue(0.0) ? tValue(0.0) : (tValue(1.0)/x);
} // safe_inverse

/** Simple shortcut for @c std::enable_shared_from_this. */
template<typename type_t>
using tObject = std::enable_shared_from_this<type_t>;

template<typename obj_t, typename obj_ptr_t>
auto shared_from_this(const obj_ptr_t& obj) {
  return std::static_pointer_cast<obj_t>(obj->shared_from_this());
} // shared_from_this

} // namespace feathers

#endif

/* TODO: Remove me. */
using int_t = feathers::ptrdiff_t;
using uint_t = feathers::size_t;
using real_t = feathers::real_t;
