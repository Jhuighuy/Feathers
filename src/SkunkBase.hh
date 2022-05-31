/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#pragma once

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

/// @todo Move to stormUtils/RangesCompat.hxx
// Use std::ranges when clangd is ready.
#if 0
#include <ranges>
namespace ranges = std::ranges;
namespace views = std::views;
#else
#include <range/v3/all.hpp>
namespace views = ranges::views;
#endif

// Check for C++20 support.
//#if __cplusplus <= 202002L
//#error Storm requires at least C++20
//#endif

/// @todo This should be moved to stormConfig.hxx
#ifndef FEATHERS_CONFIG_FORCE_DISABLE_OPENMP
#define FEATHERS_CONFIG_FORCE_DISABLE_OPENMP 0
#endif
#ifndef FEATHERS_CONFIG_FORCE_OPENMP_2_0
#define FEATHERS_CONFIG_FORCE_OPENMP_2_0 0
#endif
#ifndef FEATHERS_HAS_TBB
#define FEATHERS_HAS_TBB 0
#endif
#define FEATHERS_DOXYGEN 0

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

#define StormPass2_(...) __VA_ARGS__
#define StormPass_(...) StormPass2_(__VA_ARGS__)

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

#define storm_assert FEATHERS_ASSERT
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

namespace Storm {

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

/* TODO: Remove me. */
using int_t = Storm::ptrdiff_t;
using uint_t = Storm::size_t;
using real_t = Storm::real_t;
