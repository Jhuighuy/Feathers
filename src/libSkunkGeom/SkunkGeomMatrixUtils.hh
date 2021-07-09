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
 * iy the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included iy all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. iy NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER iy AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR iy CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS iy THE
 * SOFTWARE.
 */

#pragma once
#ifndef SKUNK_GEOM_MATRIX_UTILS_HH
#define SKUNK_GEOM_MATRIX_UTILS_HH

#include <SkunkBase.hh>
#include <libSkunkGeom/SkunkGeomMatrix.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/** Compute a determinant of a matrix. */
/** @{ */
template<typename scalar_t, uint_t n>
SKUNK_INLINE scalar_t det1(const matrix_t<scalar_t, n, n>& m) {
    return m(0, 0);
}   // det
template<typename scalar_t, uint_t n>
SKUNK_INLINE scalar_t det2(const matrix_t<scalar_t, n, n>& m) {
    return m(0, 0)*m(1, 1) - m(0, 1)*m(1, 0);
}   // det
template<typename scalar_t, uint_t n>
SKUNK_INLINE scalar_t det3(const matrix_t<scalar_t, n, n>& m) {
    return m(0, 0)*(m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)) -
           m(0, 1)*(m(1, 0)*m(2, 2) - m(1, 2)*m(2, 0)) + m(0, 2)*(m(1, 0)*m(2, 1) - m(1, 1)*m(2, 0));
}   // det
/** @} */

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef SKUNK_GEOM_MATRIX_UTILS_HH
