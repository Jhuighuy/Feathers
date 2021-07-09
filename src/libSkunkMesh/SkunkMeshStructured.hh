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
#ifndef SKUNK_MESH_STRUCTURED_HH
#define SKUNK_MESH_STRUCTURED_HH

#include <libSkunkMesh/SkunkMesh.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/** Structured mesh. */
class structured_mesh_t : public mesh_t {
public:
    /** Construct a uniform structured 1D segment mesh. */
    SKUNK_EXTERN structured_mesh_t(real_t x0, real_t x1, uint_t nx,
                                   uint_t x0_mark = 1, uint_t x1_mark = 1);
    /** Construct a uniform structured 2D rectangular mesh. */
    SKUNK_EXTERN structured_mesh_t(real_t x0, real_t x1, uint_t nx,
                                   real_t y0, real_t y1, uint_t ny,
                                   uint_t x0_mark = 1, uint_t x1_mark = 1,
                                   uint_t y0_mark = 1, uint_t y1_mark = 1,
                                   real_t skew_angle = 0.0, bool triangulate = false);
    /** Construct a uniform structured 3D hexahedral mesh. */
    SKUNK_EXTERN structured_mesh_t(real_t x0, real_t x1, uint_t nx,
                                   real_t y0, real_t y1, uint_t ny,
                                   real_t z0, real_t z1, uint_t nz,
                                   uint_t x0_mark = 1, uint_t x1_mark = 1,
                                   uint_t y0_mark = 1, uint_t y1_mark = 1,
                                   uint_t z0_mark = 1, uint_t z1_mark = 1);

    /**************************************************************************/
    /**************************************************************************/
};  // class structured_mesh_t

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef SKUNK_MESH_STRUCTURED_HH
