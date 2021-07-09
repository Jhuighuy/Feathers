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
#ifndef MHD_SPARSE_FUNCTION_HH
#define MHD_SPARSE_FUNCTION_HH

#include "SkunkBase.hh"
#include "libSkunkSparse/SkunkSparseField.hh"

template<int_t num_vars_t>
class TPiecewiseLinearFunction {
public:
    using array_t = std::array<real_t, num_vars_t>;
public:
    TScalarField<num_vars_t> u;
    TVectorField<num_vars_t> grad_u;

    TPiecewiseLinearFunction(const TScalarField<num_vars_t>& v)
        : u(v), grad_u(v.size()) {
    }

    int size() const {
        return u.size();
    }

public:
    /** Get field variables array at index. */
    /** @{ */
    array_t operator[](int_t index) const {
        array_t value{};
        for (int_t i = 0; i < num_vars_t; ++i) {
            value[i] = u[index][i];
        }
        return value;
    }
    array_t operator()(int_t index, const skunk::vec3_t& r, real_t a=1.0, real_t b=1.0) const {
        array_t value{};
        for (int_t i = 0; i < num_vars_t; ++i) {
            value[i] = a*u[index][i] + b*grad_u[index][i].inner_prod(r);
        }
        return value;
    }
    /** @} */
};

#endif  // ifndef MHD_SPARSE_FUNCTION_HH
