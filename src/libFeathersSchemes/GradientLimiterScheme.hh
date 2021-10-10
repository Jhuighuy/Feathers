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
#ifndef GRADIENT_LIMITER_SCHEME_HH_
#define GRADIENT_LIMITER_SCHEME_HH_

#include "SkunkBase.hh"
#include "libFeathersMesh/Mesh.hh"

namespace feathers {

/**
 * Barth-Jespersen (minmod)
 * slope limiter for the limiter estimation scheme.
 *
 * This is a non-differentiable limiter, so it may
 * affect convergence properties of the implicit scheme.
 */
class tMinmodSlopeLimiter {
public:
    /** Compute local slope coefficient. */
    real_t operator()(real_t du_min, real_t du_max,
                      real_t du_face, real_t eps_sqr) const;
}; // class tMinmodSlopeLimiter

/**
 * Venkatakrishnan
 * slope limiter for the limiter estimation scheme.
 *
 * This is a differentiable limiter.
 */
class tVenkatakrishnanSlopeLimiter {
public:
    /** Compute local slope coefficient. */
    real_t operator()(real_t du_min, real_t du_max,
                      real_t du_face, real_t eps_sqr) const;
}; // class tVenkatakrishnanSlopeLimiter

/**
 * Michalak Ollivier-Gooch (cubic)
 * slope limiter for the limiter estimation scheme.
 *
 * This is a differentiable limiter.
 */
class tCubicSlopeLimiter {
public:
    /** Compute local slope coefficient. */
    real_t operator()(real_t du_min, real_t du_max,
                      real_t du_face, real_t eps_sqr) const;
}; // class tCubicSlopeLimiter

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Dummy
 * second slope limiter for the limiter estimation scheme.
 *
 * This is a differentiable limiter.
 */
class tDummySecondLimiter {
public:
    /** Compute second slope coefficient. */
    real_t operator()(real_t limiter,
                      real_t du_min, real_t du_max,
                      real_t eps_sqr) const;
}; // class tDummySecondLimiter

/**
 * Michalak Ollivier-Gooch cubic
 * second slope limiter for the limiter estimation scheme.

 * This is a differentiable limiter.
 */
class tCubicSecondLimiter {
public:
    /** Compute second slope coefficient. */
    real_t operator()(real_t limiter,
                      real_t du_min, real_t du_max,
                      real_t eps_sqr) const;
}; // class tCubicSecondLimiter

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Gradient cell limiter estimation scheme:
 * computes cell-centered limiters and averages based on the cell-centered expansion.
 */
template<int_t num_vars>
class iGradientLimiterScheme : public tObject<iGradientLimiterScheme<num_vars>> {
public:
    /** Compute cell-centered gradient limiter coefficients. */
    virtual void get_cell_limiter(tScalarField<num_vars>& lim_u,
                                  const tScalarField<num_vars>& u,
                                  const tVectorField<num_vars>& grad_u) const = 0;
}; // class iGradientLimiterScheme

/**
 * Gradient cell limiter estimation scheme:
 * computes cell-centered limiters and averages based on the cell-centered expansion.
 */
template<int_t num_vars, class tSlopeLimiter, class tSecondLimiter = tDummySecondLimiter>
class tGradientLimiterScheme final : public iGradientLimiterScheme<num_vars> {
public:
    std::shared_ptr<const cMesh> m_mesh;
    tSlopeLimiter m_slope_limiter;
    tSecondLimiter m_second_limiter;

public:
    /** Initialize the limiting scheme. */
    explicit tGradientLimiterScheme(std::shared_ptr<const cMesh> mesh,
                                    const tSlopeLimiter& slope_limiter = {},
                                    const tSecondLimiter& second_limiter = {}):
        m_mesh(std::move(mesh)),
        m_slope_limiter(slope_limiter), m_second_limiter(second_limiter) {
    }

    /** Compute cell-centered gradient limiter coefficients. */
    void get_cell_limiter(tScalarField<num_vars>& lim_u,
                          const tScalarField<num_vars>& u,
                          const tVectorField<num_vars>& grad_u) const final;
}; // class iGradientLimiterScheme

template<int_t num_vars>
using tMinmodGradientLimiterScheme =
    tGradientLimiterScheme<num_vars, tMinmodSlopeLimiter>;

template<int_t num_vars>
using tVenkatakrishnanGradientLimiterScheme =
    tGradientLimiterScheme<num_vars, tVenkatakrishnanSlopeLimiter>;

template<int_t num_vars>
using tVenkatakrishnan2GradientLimiterScheme =
    tGradientLimiterScheme<num_vars, tVenkatakrishnanSlopeLimiter, tCubicSecondLimiter>;

template<int_t num_vars>
using tCubicGradientLimiterScheme =
    tGradientLimiterScheme<num_vars, tCubicSlopeLimiter>;

template<int_t num_vars>
using tCubic2GradientLimiterScheme =
    tGradientLimiterScheme<num_vars, tCubicSlopeLimiter, tCubicSecondLimiter>;

}   // namespace feathers

#include "GradientLimiterScheme.inl"

#endif // GRADIENT_LIMITER_SCHEME_HH_
