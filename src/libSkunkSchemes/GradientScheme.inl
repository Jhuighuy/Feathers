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

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/**
 * Init the gradient scheme.
 */
template<int_t num_vars>
void TLeastSquaresGradientScheme<num_vars>::init_gradients_() {
    /* Compute the least-squares
     * problem matrices for the interior cells. */
    for_each_interior_cell(*m_mesh, [&](CellIter cell) {
        mat3_t& mat = m_inverse_matrices[cell][0];
        cell.for_each_face_cells([&](CellIter cell_inner, CellIter cell_outer) {
            const vec3_t dr =
                cell_outer->get_center_position() - cell_inner->get_center_position();
            mat += dr.outer_prod(dr);
        });
    });

    /* Compute the least squares problem
     * right hand statements for the boundary cells.
     * Use the same stencil as for the interior cell, but centered to a boundary cell. */
    for_each_boundary_face_cells(*m_mesh, [&](CellIter cell_inner, CellIter cell_outer) {
        mat3_t& mat = m_inverse_matrices[cell_outer][0];
        const vec3_t dr =
            cell_outer->get_center_position() - cell_inner->get_center_position();
        mat += dr.outer_prod(dr);
        cell_inner.for_each_face_cells([&](CellIter cell_inner_inner,
                                           CellIter cell_inner_outer) {
            if (cell_inner_inner == cell_inner) {
                const vec3_t dr_inner =
                    cell_inner_outer->get_center_position() - cell_inner->get_center_position();
                mat += dr_inner.outer_prod(dr_inner);
            } else if (cell_inner_outer == cell_inner) {
                const vec3_t dr_inner =
                    cell_inner_inner->get_center_position() - cell_inner->get_center_position();
                mat += dr_inner.outer_prod(dr_inner);
            }
        });
    });

    /* Compute the inverse of the least squares problem matrices.
     * ( Matrix is stabilized by a small number, added to the diagonal. ) */
    for_each_cell(*m_mesh, [&](CellIter cell) {
        static const mat3_t eps{1e-14, 0.0, 0.0,
                                0.0, 1e-14, 0.0,
                                0.0, 0.0, 1e-14};
        mat3_t& mat = m_inverse_matrices[cell][0];
        mat = mat3_t::inv(mat + eps);
    });
}   // TLeastSquaresGradientScheme::init_gradients_

/**
 * Compute cell-centered gradients using the Weighted Least-Squares, cell-based version.
 */
template<int_t num_vars>
template</*template<int_t>*/ class TOutField,
         /*template<int_t>*/ class TInField>
void TLeastSquaresGradientScheme<num_vars>::get_gradients_(TOutField/*<num_vars>*/& grad_u,
                                                                const TInField/*<num_vars>*/& u) const {
    /* Compute the least-squares
     * problem right hand statements for the interior cells. */
    for_each_interior_cell(*m_mesh, [&](CellIter cell) {
        cell.for_each_face_cells([&](CellIter cell_inner, CellIter cell_outer) {
          const vec3_t dr =
              cell_outer->get_center_position() - cell_inner->get_center_position();
            for (int_t i = 0; i < num_vars; ++i) {
                grad_u[cell][i] += (u[cell_outer][i] - u[cell_inner][i])*dr;
            }
        });
    });

    /* Compute the least squares problem
     * right hand statements for the boundary cells.
     * Use the same stencil as for the interior cell, but centered to a boundary cell. */
    for_each_boundary_face_cells(*m_mesh, [&](CellIter cell_inner, CellIter cell_outer) {
      const vec3_t dr =
          cell_outer->get_center_position() - cell_inner->get_center_position();
        for (int_t i = 0; i < num_vars; ++i) {
            grad_u[cell_outer][i] += (u[cell_outer][i] - u[cell_inner][i])*dr;
        }
        cell_inner.for_each_face_cells([&](CellIter cell_inner_inner,
                                           CellIter cell_inner_outer) {
            if (cell_inner_inner == cell_inner) {
                const vec3_t dr_inner =
                    cell_inner_outer->get_center_position() - cell_inner->get_center_position();
                for (int_t i = 0; i < num_vars; ++i) {
                    grad_u[cell_outer][i] += (u[cell_inner_outer][i] - u[cell_inner][i])*dr_inner;
                }
            } else if (cell_inner_outer == cell_inner) {
                const vec3_t dr_inner =
                    cell_inner_inner->get_center_position() - cell_inner->get_center_position();
                for (int_t i = 0; i < num_vars; ++i) {
                    grad_u[cell_outer][i] += (u[cell_inner_inner][i] - u[cell_inner][i])*dr_inner;
                }
            }
        });
    });

    /* Solve the least-squares problem. */
    for_each_cell(*m_mesh, [&](CellIter cell) {
        for (int_t i = 0; i < num_vars; ++i) {
            const mat3_t& mat = m_inverse_matrices[cell][0];
            grad_u[cell][i] = mat.prod(grad_u[cell][i]);
        }
    });
}   // TLeastSquaresGradientScheme::get_gradients_

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
