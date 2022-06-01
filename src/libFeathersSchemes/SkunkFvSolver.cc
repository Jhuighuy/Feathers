// ************************************************************************************
// // Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver. Copyright(C)
// Butakov Oleg and Co. 2019.
// ************************************************************************************
// //

#include "SkunkFvSolver.hh"

// ************************************************************************************
// //
// ************************************************************************************
// //
// ************************************************************************************
// //

template<typename MhdPhysicsT>
MhdFvSolverT<MhdPhysicsT>::MhdFvSolverT(std::shared_ptr<const cMesh> mesh) :
    m_mesh(mesh), m_conv(new Storm::cUpwind2ConvectionScheme(mesh)) {
  m_bcs[1] = std::make_shared<MhdFvBcFarFieldT<MhdPhysicsT>>();
  m_bcs[2] = std::make_shared<MhdFvBcSlipT<MhdPhysicsT>>();
}

// ------------------------------------------------------------------------------------
// //
// ------------------------------------------------------------------------------------
// //

/**
 * @brief Compute spacial discretization.
 */
template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_func(Storm::tScalarField& u,
                                          Storm::tScalarField& u_out) const {
  using namespace Storm;

  /*
   * Clear fields and apply boundary conditions.
   */

  ForEach(CellViews(*m_mesh), [&](CellView cell) { u_out[cell].fill(0.0); });
  for (size_t mark = 1; mark < m_mesh->NumFaceMarks(); ++mark) {
    const auto& bc = m_bcs.at(mark);
    ForEach(FaceViews(*m_mesh, FaceMark(mark)), [&](FaceView face) {
      bc->get_ghost_state(face.Normal(), face.InnerCell().Center(),
                          face.OuterCell().Center(), u[face.InnerCell()].data(),
                          u[face.OuterCell()].data());
    });
  }

  m_conv->get_cell_convection(5, u_out, u);

} // MhdFvSolverT::calc_func

template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_step(real_t& dt, Storm::tScalarField& u,
                                          Storm::tScalarField& u_hat) const {
  /*
   * Compute.
   */
  calc_func(u, u_hat);
  ForEach(IntCellViews(*m_mesh), [&](CellView cell) {
    for (uint_t i = 0; i < num_vars; ++i) {
      u_hat[cell][i] = u[cell][i] - dt * u_hat[cell][i];
    }
  });
} // MhdFvSolverT::calc_step

// ************************************************************************************
// //
// ************************************************************************************
// //
// ************************************************************************************
// //

template class MhdFvSolverT<tGasPhysics>;
