// Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver.
// Copyright(C) Butakov Oleg and Co. 2019.

#pragma once

#include "SkunkBase.hh"
#include <cstring>

// ************************************************************************************
// //
// ************************************************************************************
// //
// ************************************************************************************
// //

static const real_t Gamma = 1.4;
static const real_t Gamma1 = Gamma - 1.0;

using namespace Storm;

class MhdHydroVars {
public:

  real_t rho = 0.0; /**< Fluid density, ๐. */
  real_t p = 0.0;   /**< Fluid pressure, ๐. */
  vec3_t vel = {};  /**< Fluid velocity, ๐. */
  real_t vel_n = 0.0; /**< Fluid velocity normal component, ๐โ = ๐โ๐. */
  real_t eps = 0.0; /**< Fluid internal energy, ๐. */
  real_t nrg = 0.0; /**< Fluid specific total energy, ๐ธ = ยฝ๐ยฒ + ๐. */
  real_t ent = 0.0; /**< Fluid specific enthalpy, ๐ป = ๐ธ + ๐/๐. */
  real_t c_snd = 0.0; /**< Fluid sound speed, ๐ = (โ๐/โ๐)ยนแยฒ. */
  const real_t* rest_prim = nullptr; /**< Additional advected scalars, ๐แตข. */
  const real_t* rest_cons =
      nullptr; /**< Additional advected scalars, ๐ขแตข = ๐๐แตข. */

public:

  explicit MhdHydroVars() = default;
  explicit MhdHydroVars(const vec3_t& n, const real_t* q_cons,
                        const real_t* q_prim = nullptr);

public:

  /** make primitive variables, ๐ธ = (๐,๐,๐,๐แตข,โฆ)แต. */
  real_t* make_prim(uint_t num_vars, real_t* prim) const {
    *reinterpret_cast<std::array<real_t, 5>*>(prim) = {rho, p, vel.x, vel.y,
                                                       vel.z};
    return prim;
  }

  /** make conserved variables, ๐ผ = (๐,๐๐ธ,๐๐,๐๐แตข,โฆ)แต. */
  real_t* make_cons(uint_t num_vars, real_t* cons) const {
    *reinterpret_cast<std::array<real_t, 5>*>(cons) = {
        rho, rho * nrg, rho * vel.x, rho * vel.y, rho * vel.z};
    for (uint_t i = 5; i < num_vars; ++i) {
      if (rest_cons != nullptr) {
        cons[i] = rest_cons[i - 5];
      } else if (rest_prim != nullptr) {
        cons[i] = rho * rest_prim[i - 5];
      }
    }
    return cons;
  }

  /** make flux variables, ๐ญโ = (๐๐โ,๐๐ป๐โ,๐๐๐โ,๐๐แตข๐โ,โฆ)แต. */
  real_t* make_flux(uint_t num_vars, const vec3_t& n, real_t* flux) const {
    *reinterpret_cast<std::array<real_t, 5>*>(flux) = {
        rho * vel_n, rho * vel_n * ent, rho * vel_n * vel.x + p * n.x,
        rho * vel_n * vel.y + p * n.y, rho * vel_n * vel.z + p * n.z};
    for (uint_t i = 5; i < num_vars; ++i) {
      if (rest_cons != nullptr) {
        flux[i] = vel_n * rest_cons[i - 5];
      } else if (rest_prim != nullptr) {
        flux[i] = rho * vel_n * rest_prim[i - 5];
      }
    }
    return flux;
  }
};

inline MhdHydroVars::MhdHydroVars(const vec3_t& n, const real_t* q_cons,
                                  const real_t* q_prim)
    : MhdHydroVars() {
  if (q_cons != nullptr) {
    rho = q_cons[0];
    nrg = q_cons[1] / rho;
    vel.x = q_cons[2] / rho;
    vel.y = q_cons[3] / rho;
    vel.z = q_cons[4] / rho;
    vel_n = glm::dot(vel, n);
    eps = nrg - 0.5 * glm::dot(vel, vel);
    p = Gamma1 * rho * eps;
    ent = nrg + p / rho;
  } else if (q_prim != nullptr) {
    rho = q_prim[0];
    p = q_prim[1];
    vel.x = q_prim[2];
    vel.y = q_prim[3];
    vel.z = q_prim[4];
    vel_n = glm::dot(vel, n);
    eps = p / rho / Gamma1;
    nrg = eps + 0.5 * glm::dot(vel, vel);
    ent = nrg + p / rho;
  }
  c_snd = std::sqrt(Gamma * p / rho);
}

typedef class MhdHydroVars MhdFluidVarsIdealGas;

class tGasPhysics {
public:

  static constexpr int_t num_vars = 5;
  typedef MhdFluidVarsIdealGas MhdFluidStateT;
  typedef MhdFluidVarsIdealGas tFluidState;
}; // class tGasPhysics
