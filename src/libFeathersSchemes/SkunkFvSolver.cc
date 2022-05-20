// ************************************************************************************ //
// Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver.
// Copyright(C) Butakov Oleg and Co. 2019.
// ************************************************************************************ //

#include "SkunkFvSolver.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template<typename MhdPhysicsT>
MhdFvSolverT<MhdPhysicsT>::MhdFvSolverT(std::shared_ptr<const cMesh> mesh)
    : m_mesh(mesh),
      m_conv(new feathers::cUpwind2ConvectionScheme(mesh)) {
    m_bcs[1] = std::make_shared<MhdFvBcFarFieldT<MhdPhysicsT>>();
    m_bcs[2] = std::make_shared<MhdFvBcSlipT<MhdPhysicsT>>();
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * @brief Compute spacial discretization.
 */
template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_func(feathers::tScalarField& u,
                                          feathers::tScalarField& u_out) const {

  using namespace feathers;

  /*
   * Clear fields and apply boundary conditions.
   */

  ForEach(CellRefs(*m_mesh), [&](CellRef cell) {
    u_out[cell].fill(0.0);
  });
  for (size_t mark = 1; mark < m_mesh->NumFaceMarks(); ++mark) {
      const auto& bc = m_bcs.at(mark);
    ForEach(FaceRefs(*m_mesh, FaceMark(mark)), [&](FaceRef face) {
      bc->get_ghost_state(face.Normal(),
                          face.InnerCell().CenterPos(),
                          face.OuterCell().CenterPos(),
                          u[face.InnerCell()].data(),
                          u[face.OuterCell()].data());
    });
  }

  m_conv->get_cell_convection(5, u_out, u);

}   // MhdFvSolverT::calc_func

template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_step(real_t& dt,
                                          feathers::tScalarField& u,
                                          feathers::tScalarField& u_hat) const {
  /*
   * Compute.
   */
  calc_func(u, u_hat);
  ForEach(InteriorCellRefs(*m_mesh), [&](CellRef cell) {
    for (uint_t i = 0; i < num_vars; ++i) {
      u_hat[cell][i] = u[cell][i] - dt * u_hat[cell][i];
    }
  });
}   // MhdFvSolverT::calc_step

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template class MhdFvSolverT<tGasPhysics>;
