// Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver.
// Copyright(C) Butakov Oleg and Co. 2019.

#pragma once

#include "SkunkBase.hh"
#include <cstring>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

static const real_t Gamma = 1.4;
static const real_t Gamma1 = Gamma - 1.0;

using namespace feathers;

class MhdHydroVars {
public:
    real_t rho   = 0.0; /**< Fluid density, 𝜌. */
    real_t p     = 0.0; /**< Fluid pressure, 𝑝. */
    vec3_t vel   = { }; /**< Fluid velocity, 𝒗. */
    real_t vel_n = 0.0; /**< Fluid velocity normal component, 𝒗ₙ = 𝒗⋅𝒏. */
    real_t kin   = 0.0; /**< Fluid specific kinetic energy, 𝐾 = ½𝒗². */
    real_t eps   = 0.0; /**< Fluid internal energy, 𝜀. */
    real_t nrg   = 0.0; /**< Fluid specific total energy, 𝐸 = 𝐾 + 𝜀. */
    real_t ent   = 0.0; /**< Fluid specific enthalpy, 𝐻 = 𝐸 + 𝑝/𝜌. */
    real_t c_snd = 0.0; /**< Fluid sound speed, 𝑐 = (∂𝑝/∂𝜌)¹ᐟ². */

    std::array<real_t, 5> prim = {}; /**< Primitive variables, 𝑷 = (𝜌,𝑝,𝒗)ᵀ. */
    std::array<real_t, 5> cons = {}; /**< Conserved variables, 𝑼 = (𝜌,𝜌𝐸,𝜌𝒗)ᵀ. */
    std::array<real_t, 5> flux = {}; /**< Flux variables, 𝑭ₙ = (𝜌𝒗ₙ,𝜌𝐻𝒗ₙ,𝜌𝒗𝒗ₙ)ᵀ. */

public:
    explicit MhdHydroVars() = default;
    explicit MhdHydroVars(const vec3_t& n,
                          const real_t* q_cons,
                          const real_t* q_prim = nullptr);

public:
    void make_cons() {
        cons = { rho, rho*nrg, rho*vel.x, rho*vel.y, rho*vel.z };
    }
};

inline
MhdHydroVars::MhdHydroVars(const vec3_t& n,
                           const real_t* q_cons, const real_t* q_prim) : MhdHydroVars() {
    if (q_cons != nullptr) {
        rho = q_cons[0];
        nrg = q_cons[1]/rho;
        vel.x = q_cons[2]/rho;
        vel.y = q_cons[3]/rho;
        vel.z = q_cons[4]/rho;
        vel_n = glm::dot(vel, n);
        kin = 0.5*glm::dot(vel, vel);
        eps = nrg - kin;
        p   = Gamma1*rho*eps;
        ent = nrg + p/rho;
    } else if (q_prim != nullptr) {
        rho = q_prim[0];
        p   = q_prim[1];
        vel.x = q_prim[2];
        vel.y = q_prim[3];
        vel.z = q_prim[4];
        vel_n = glm::dot(vel, n);
        kin = 0.5*glm::dot(vel, vel);
        eps = p/rho/Gamma1;
        nrg = eps + kin;
        ent = nrg + p/rho;
    }

    c_snd = std::sqrt(Gamma*p/rho);
    prim = { rho, p, vel.x, vel.y, vel.z };
    cons = { rho, rho*nrg, rho*vel.x, rho*vel.y, rho*vel.z };
    flux = { rho*vel_n, rho*vel_n* ent, rho*vel_n*vel.x + p*n.x,
             rho*vel_n*vel.y + p*n.y, rho*vel_n*vel.z + p*n.z };
}

typedef class MhdHydroVars MhdFluidVarsIdealGas;

class tGasPhysics {
public:
    static constexpr int_t num_vars = 5;
    typedef MhdFluidVarsIdealGas MhdFluidStateT;
}; // class tGasPhysics
