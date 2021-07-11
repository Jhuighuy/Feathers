// ************************************************************************************ //
// Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver.
// Copyright(C) Butakov Oleg and Co. 2019.
// ************************************************************************************ //

#include "SkunkFvSolver.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template<typename MhdPhysicsT>
MhdFvSolverT<MhdPhysicsT>::MhdFvSolverT(std::shared_ptr<cMesh> mesh)
    : m_mesh(mesh),
      m_conv(new feathers::tUpwindConvectionScheme<num_vars>(mesh)) {
    m_bcs[1] = std::make_shared<MhdFvBcFarFieldT<MhdPhysicsT>>();
    m_bcs[2] = std::make_shared<MhdFvBcSlipT<MhdPhysicsT>>();
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * @brief Compute spacial discretization.
 */
template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_func(tScalarField<num_vars>& u,
                                          tScalarField<num_vars>& u_out) const {
    /*
     * Clear fields and apply boundary conditions.
     */
#pragma omp parallel for
    for (uint_t cell_ind = 0; cell_ind < m_mesh->num_cells(); ++cell_ind) {
        u_out[cell_ind].fill(0.0);
    }
    for (uint_t mark = 1; mark < m_mesh->num_face_marks(); ++mark) {
#pragma omp parallel for
        for (uint_t face_ind = m_mesh->begin_face(mark); face_ind != m_mesh->end_face(mark); ++face_ind) {
            const auto& face = m_mesh->get_face(face_ind);
            m_bcs.at(mark)->get_ghost_state(m_mesh->get_face_normal(face_ind),
                                            m_mesh->get_cell_center_position(face.get_inner_cell()),
                                            m_mesh->get_cell_center_position(face.get_outer_cell()),
                                            u[face.get_inner_cell()],
                                            u[face.get_outer_cell()]);
        }
    }
    m_conv->get_cell_convection(u_out, u);
}   // MhdFvSolverT::calc_func

template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_step(real_t& dt,
                                          tScalarField<num_vars>& uc,
                                          tScalarField<num_vars>& up) const {
    /*
     * Compute.
     */
    calc_func(uc, up);
#pragma omp parallel for
    for (uint_t cell_ind = m_mesh->begin_cell(0); cell_ind != m_mesh->end_cell(0); ++cell_ind) {
        for (uint_t i = 0; i < num_vars; ++i) {
            up[cell_ind][i] = uc[cell_ind][i] - up[cell_ind][i]*dt;
        }
    }
}   // MhdFvSolverT::calc_step

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template class MhdFvSolverT<MhdPhysicsIdealGas>;
