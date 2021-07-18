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

#include "ConvectionScheme.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Compute the first-order upwind convection.
 */
template<int_t num_vars>
void tUpwindConvectionScheme<num_vars>::get_cell_convection(tScalarField<num_vars>& div_f,
                                                            const tScalarField<num_vars>& u) const {
    /* Compute the first order numerical fluxes. */
    tScalarField<num_vars> flux_u(m_mesh->num_faces());
    for_each_face(*m_mesh, [&](tFaceIter face) {
        const tCellIter cell_outer = face.get_outer_cell();
        const tCellIter cell_inner = face.get_inner_cell();

        std::array<real_t, num_vars>& flux = flux_u[face];
        flux.fill(0.0), m_flux->get_numerical_flux(face.get_normal(),
                                                   u[cell_outer], u[cell_inner], flux);
    });

    /* Compute the first order convection. */
    for_each_interior_cell(*m_mesh, [&](tCellIter cell) {
        std::array<real_t, num_vars>& conv = div_f[cell];
        conv.fill(0.0), cell.for_each_face([&](tFaceIter face) {
            const tCellIter cell_outer = face.get_outer_cell();
            const tCellIter cell_inner = face.get_inner_cell();
            const real_t ds = face.get_area();
            if (cell_outer == cell) {
                for (int_t i = 0; i < num_vars; ++i) {
                    conv[i] -= flux_u[face][i] * ds;
                }
            } else if (cell_inner == cell) {
                for (int_t i = 0; i < num_vars; ++i) {
                    conv[i] += flux_u[face][i] * ds;
                }
            }
        });
        const real_t inv_dv = 1.0/cell.get_volume();
        for (int_t i = 0; i < num_vars; ++i) {
            conv[i] *= inv_dv;
        }
    });
}   // tUpwindConvectionScheme<num_vars>::get_cell_convection

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/**
 * Compute the second-order upwind convection.
 */
template<int_t num_vars>
void tUpwind2ConvectionScheme<num_vars>::get_cell_convection(tScalarField<num_vars>& div_f,
                                                             const tScalarField<num_vars>& u) const {
    /* Compute the second order limited gradients. */
    tPiecewiseLinearFunction<num_vars> u_rec(u);
    m_gradient_scheme->get_gradients(u_rec.grad_u, u);

    tScalarField<num_vars> lim_u(m_mesh->num_cells());
    m_gradient_limiter_scheme->get_cell_limiter(lim_u, u_rec);

    for_each_cell(*m_mesh, [&](tCellIter cell) {
        for (int_t i = 0; i < num_vars; ++i) {
            u_rec.grad_u[cell][i] *= lim_u[cell][i];
        }
    });

    /* Compute the second order numerical fluxes:
     * integrate the numerical flux over the face nodes. */
    tScalarField<num_vars> flux_f(m_mesh->num_faces());
    for_each_face(*m_mesh, [&](tFaceIter face) {
        const tCellIter cell_outer = face.get_outer_cell();
        const tCellIter cell_inner = face.get_inner_cell();
        const vec3_t dr_outer =
            face.get_center_coords() - cell_outer.get_center_coords();
        const vec3_t dr_inner =
            face.get_center_coords() - cell_inner.get_center_coords();
        const auto u_outer = u_rec(cell_outer, dr_outer);
        const auto u_inner = u_rec(cell_inner, dr_inner);

        std::array<real_t, num_vars>& flux = flux_f[face];
        flux.fill(0.0), m_flux->get_numerical_flux(face.get_normal(), u_outer, u_inner, flux);
    });

    /* Compute the second order convection. */
    for_each_interior_cell(*m_mesh, [&](tCellIter cell) {
        div_f[cell].fill(0.0);
        cell.for_each_face([&](tFaceIter face) {
            const tCellIter cell_outer = face.get_outer_cell();
            const tCellIter cell_inner = face.get_inner_cell();
            const real_t ds = face.get_area();
            if (cell_outer == cell) {
                for (int_t i = 0; i < num_vars; ++i) {
                    div_f[cell][i] -= flux_f[face][i] * ds;
                }
            } else if (cell_inner == cell) {
                for (int_t i = 0; i < num_vars; ++i) {
                    div_f[cell][i] += flux_f[face][i] * ds;
                }
            }
        });
        const real_t inv_dv = 1.0/cell.get_volume();
        for (int_t i = 0; i < num_vars; ++i) {
            div_f[cell][i] *= inv_dv;
        }
    });
}   // tUpwind2ConvectionScheme::get_cell_convection

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template class feathers::tUpwindConvectionScheme<5>;
template class feathers::tUpwind2ConvectionScheme<5>;
