// ************************************************************************************ //
// Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver.
// Copyright(C) Butakov Oleg and Co. 2019.
// ************************************************************************************ //

#pragma once

#include "ConvectionScheme.hh"
#include "FluxScheme.hh"
#include "GradientLimiterScheme.hh"
#include "SkunkBase.hh"
#include "SkunkFvBC.hh"

#include <map>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/**
 * @brief A Finite Volume solver.
 */
template<typename MhdPhysicsT>
class MhdFvSolverT :
    public std::enable_shared_from_this<MhdFvSolverT<MhdPhysicsT>> {
public:
    using MhdFluidStateT = typename MhdPhysicsT::MhdFluidStateT;
    static constexpr int_t num_vars = MhdPhysicsT::num_vars;

private:
    std::shared_ptr<UMesh> m_mesh;
    std::shared_ptr<feathers::TConvectionScheme<num_vars>> m_conv;
    std::map<int_t, std::shared_ptr<MhdFvBcPT<MhdPhysicsT>>> m_bcs;

public:
    explicit MhdFvSolverT(std::shared_ptr<UMesh> mesh);

public:
    void calc_func(TScalarField<num_vars>& u,
                   TScalarField<num_vars>& up) const;
    void calc_step(real_t& dt,
                   TScalarField<num_vars>& u,
                   TScalarField<num_vars>& up) const;
};  // class MhdFvSolverT

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
