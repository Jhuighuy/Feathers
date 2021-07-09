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
#ifndef GRADIENT_SCHEME_HH_
#define GRADIENT_SCHEME_HH_

#include "SkunkBase.hh"
#include "libSkunkMesh/SkunkMesh.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/** Abstract cell-centered gradient scheme. */
template<int_t num_vars>
class IGradientScheme : public TObject<IGradientScheme<num_vars>> {
public:
    /** Compute cell-centered gradients. */
    /** @{ */
    virtual void get_gradients(TVectorField<num_vars>& grad_u,
                               const TScalarField<num_vars>& u) const = 0;
    /** @} */
};  // class IGradientScheme

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/**
 * Weighted Least-Squares gradient estimation scheme, cell-based:
 * computes cell-centered gradients based on the cell-centered values.
 *
 * This gradient scheme is a second-order scheme for any meshes.
 * Also, this gradient scheme is by far the fastest one.
 */
template<int_t num_vars>
class TLeastSquaresGradientScheme final : public IGradientScheme<num_vars> {
private:
    std::shared_ptr<mesh_t> m_mesh;
    TMatrixField<1> m_inverse_matrices;

public:
    /** Initialize the gradient scheme. */
    explicit TLeastSquaresGradientScheme(std::shared_ptr<UMesh> mesh):
        m_mesh(std::move(mesh)),
        m_inverse_matrices(m_mesh->num_cells()) {
        init_gradients_();
    }

    /** Compute cell-centered gradients. */
    /** @{ */
    void get_gradients(TVectorField<num_vars>& grad_u,
                       const TScalarField<num_vars>& u) const final {
        get_gradients_(grad_u, u);
    }
    /** @} */

private:
    /** Compute cell-centered gradients. */
    /** @{ */
    void init_gradients_();
    template</*template<int_t>*/ class TOutField,
             /*template<int_t>*/ class TInField>
    void get_gradients_(TOutField/*<num_vars>*/& grad_u,
                        const TInField/*<num_vars>*/& u) const;
    /** @} */
};  // class TLeastSquaresGradientScheme

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#include "GradientScheme.inl"

#endif // GRADIENT_SCHEME_HH_
