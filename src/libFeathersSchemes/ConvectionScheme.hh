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
#ifndef CONVECTION_SCHEME_HH_
#define CONVECTION_SCHEME_HH_

#include "SkunkBase.hh"
#include "FluxScheme.hh"
#include "GradientLimiterScheme.hh"
#include "GradientScheme.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Abstract convection scheme.
 */
template<int_t num_vars>
class iConvectionScheme : public tObject<iConvectionScheme<num_vars>> {
public:
    virtual void get_cell_convection(tScalarField<num_vars>& conv_u,
                                     const tScalarField<num_vars>& u) const = 0;
};  // class iConvectionScheme

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Piecewise-constant upwind convection scheme.
 * This is a first-order scheme.
 */
template<int_t num_vars>
class tUpwindConvectionScheme final : public iConvectionScheme<num_vars> {
public:
    std::shared_ptr<const uMesh> m_mesh;
    std::shared_ptr<iFluxScheme<num_vars>> m_flux;

public:
    explicit tUpwindConvectionScheme(std::shared_ptr<const uMesh> mesh):
        m_mesh(std::move(mesh)),
        m_flux(new tHLLCFluxScheme<MhdPhysicsIdealGas>()) {
    }

    /** Compute the first-order upwind convection. */
    void get_cell_convection(tScalarField<num_vars>& div_f,
                             const tScalarField<num_vars>& u) const final;
};  // class tUpwindConvectionScheme

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Piecewise-linear upwind convection scheme.
 * This is a second-order scheme.
 */
template<int_t num_vars>
class tUpwind2ConvectionScheme final : public iConvectionScheme<num_vars> {
public:
    std::shared_ptr<const uMesh> m_mesh;
    std::shared_ptr<iFluxScheme<num_vars>> m_flux;
    std::shared_ptr<iGradientScheme<num_vars>> m_gradient_scheme;
    std::shared_ptr<iGradientLimiterScheme<num_vars>> m_gradient_limiter_scheme;

public:
    explicit tUpwind2ConvectionScheme(std::shared_ptr<const uMesh> mesh):
        m_mesh(std::move(mesh)),
        m_flux(new tLaxFriedrichsFluxScheme<MhdPhysicsIdealGas>()),
        m_gradient_scheme(new tLeastSquaresGradientScheme<num_vars>(m_mesh)),
        m_gradient_limiter_scheme(new tCubic2GradientLimiterScheme<num_vars>(m_mesh)) {
    }

    /** Compute the second-order upwind convection. */
    void get_cell_convection(tScalarField<num_vars>& div_f,
                             const tScalarField<num_vars>& u) const final;
};  // class tUpwindConvectionScheme

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // CONVECTION_SCHEME_HH_