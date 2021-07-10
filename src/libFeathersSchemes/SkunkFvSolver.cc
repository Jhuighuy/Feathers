// ************************************************************************************ //
// Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver.
// Copyright(C) Butakov Oleg and Co. 2019.
// ************************************************************************************ //

#include "SkunkFvSolver.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template<typename MhdPhysicsT>
MhdFvSolverT<MhdPhysicsT>::MhdFvSolverT(std::shared_ptr<UMesh> mesh)
    : m_mesh(mesh),
      m_conv(new feathers::tUpwind2ConvectionScheme<num_vars>(mesh)) {
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
                                          tScalarField<num_vars>& up) const {
    /*
     * Clear fields and apply boundary conditions.
     */
#pragma omp parallel for
    for (uint_t cell_ind = 0; cell_ind < m_mesh->num_cells(); ++cell_ind) {
        up[cell_ind].fill(0.0);
    }
    for (uint_t mark = 1; mark < m_mesh->num_face_marks(); ++mark) {
#pragma omp parallel for
        for (uint_t face_mind = 0; face_mind < m_mesh->num_marked_faces(mark); ++face_mind) {
            const uint_t face_ind = m_mesh->get_marked_face_index(face_mind, mark);
            const auto& face = m_mesh->get_face(face_ind);
            m_bcs.at(mark)->get_ghost_state(face.get_normal(),
                                            m_mesh->get_cell_center_position(face.get_inner_cell()),
                                            m_mesh->get_cell_center_position(face.get_outer_cell()),
                                            u[face.get_inner_cell()],
                                            u[face.get_outer_cell()]);
        }
    }

    m_conv->get_cell_convection(up, u);

#if 0
    const real_t r1 = 1.3, r2 = 3.5;
    const auto rho_bar = 0.55;
    const auto p_bar = 1.3;
    const auto r_bar = r2;
    const auto k = p_bar*std::pow(rho_bar, -Gamma);
    const auto g_bar = 1.04;
    const auto G_bar = 0.5*g_bar*r_bar;
    const auto c_bar = Gamma/Gamma1*p_bar/rho_bar + G_bar;
    const auto omega_bar = 2.306e-2*50;
    foreach_cell_interior_index(MhdMeshCellIndex, cell, m_mesh) {
        const auto x = cell->get_center().x;
        const auto y = cell->get_center().y;
        const auto z = cell->get_center().z;

#if 0
        const auto r = std::sqrt(x*x + y*y + z*z);
        const auto phi = std::atan2(y, x);
        const auto psy = std::atan2(z, std::hypot(x, y));

        const mhd_mat3_t R{
            {  cos(phi)*cos(psy), sin(phi)*cos(psy), +sin(psy) },
            {  cos(phi)*sin(psy), sin(phi)*sin(psy), -cos(psy) },
            { -sin(phi),          cos(phi),           0.0      },
        };

        const auto rho = u[cell][0];
        auto v = mhd_vec3_t{u[cell][2], u[cell][3], u[cell][4]}/rho;
        v = R*v;
        const auto v_r = v.x;
        const auto v_psy = v.y;
        const auto v_phi = v.z;

        mhd_vec3_t F{
            +rho*cos(psy)*(omega_bar*omega_bar*r*cos(psy) + 2*omega_bar*v_phi),
            -rho*sin(psy)*(omega_bar*omega_bar*r*cos(psy) + 2*omega_bar*v_phi),
            rho*omega_bar*2.0*(v_psy*sin(psy) - v_r*cos(psy)),
        };
        mhd_vec3_t G{
            rho*g_bar*r/r_bar,
            0.0,
            0.0,
        };
#endif

        const auto rho = u[cell][0];
        auto v = mhd_vec3_t{u[cell][2], u[cell][3], u[cell][4]}/rho;
        mhd_vec3_t R = {x, y, z};
        const auto r = std::sqrt(x*x + y*y + z*z);
        mhd_vec3_t g = -g_bar*R/r_bar;

        mhd_vec3_t w{0, 0, 1.0};
        auto Fc = (-mhd_vec3_t::cross(w, mhd_vec3_t::cross(w, R)) - 2.0*mhd_vec3_t::cross(w, v));

        up[cell][1] -= 0*rho*v*(g + Fc);
        ((mhd_vec3_t&)up[cell][2]) -= 0*rho*(g + Fc);
    }
#endif
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
    for (uint_t cell_mind = 0; cell_mind < m_mesh->num_marked_cells(0); ++cell_mind) {
        const uint_t cell_ind = m_mesh->get_marked_cell_index(cell_mind, 0);
        for (uint_t i = 0; i < num_vars; ++i) {
            up[cell_ind][i] = uc[cell_ind][i] - up[cell_ind][i]*dt;
        }
    }
}   // MhdFvSolverT::calc_step

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template class MhdFvSolverT<MhdPhysicsIdealGas>;
