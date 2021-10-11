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
#ifndef FIELD_HH_
#define FIELD_HH_

#include "SkunkBase.hh"

template<int num_vars = 1>
using tScalarField = std::vector<std::array<feathers::real_t, num_vars>>;
template<int num_vars = 1>
using tVectorField = std::vector<std::array<feathers::vec3_t, num_vars>>;
template<int num_vars = 1>
using tMatrixField = std::vector<std::array<feathers::mat3_t, num_vars>>;

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {



} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // FIELD_HH_