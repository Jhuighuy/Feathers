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
    real_t rho;
    real_t nrg, eps, kin;
    real_t ent, p;
    vec3_t vel;
    real_t vel_n, V2;
    real_t c_snd, c2snd;
    std::array<real_t, 5> prim;
    std::array<real_t, 5> cons;
    std::array<real_t, 5> flux;

public:
    explicit MhdHydroVars() {
        std::memset((void*)this, 0, sizeof(*this));
    }
    explicit MhdHydroVars(const vec3_t& n,
                          const real_t* q_cons,
                          const real_t* q_prim = nullptr);
public:
    void make_cons() {
        //real_t pr[]{rho, p, vel.x, vel.y, vel.z};
        //cons = MhdHydroVars({}, nullptr, pr).cons;
#if 1
        cons = {
                rho,
                rho*nrg,
            rho * vel.x, rho * vel.y, rho * vel.z,
        };
#endif
    }
};

inline
MhdHydroVars::MhdHydroVars(const vec3_t& n,
                           const real_t* q_cons, const real_t* q_prim) : MhdHydroVars() {
    if (q_cons != nullptr) {
        rho = q_cons[0];
        if (rho == 0.0) return;
        nrg = q_cons[1]/rho;
        vel.x = q_cons[2] / rho;
        vel.y = q_cons[3] / rho;
        vel.z = q_cons[4] / rho;
        vel_n  = glm::dot(vel, n);
        V2  = glm::dot(vel, vel);
        kin = 0.5*V2;
        eps = nrg - kin;
        p   = Gamma1*rho*eps;
        ent = nrg + p/rho;
    } else if (q_prim != nullptr) {
        rho = q_prim[0];
        p   = q_prim[1];
        vel.x = q_prim[2];
        vel.y = q_prim[3];
        vel.z = q_prim[4];
        vel_n  = glm::dot(vel, n);
        V2  = glm::dot(vel, vel);
        kin = 0.5*V2;
        if (rho == 0.0) return;
        eps = p/rho/Gamma1;
        nrg = eps + kin;
        ent = nrg + p/rho;
    }
    c2snd = Gamma*p/rho;
    c_snd = std::sqrt(c2snd);
    if (std::isnan(c_snd)) {
        abort();
    }
    prim = {
        rho, p, vel.x, vel.y, vel.z,
    };
    cons = {
        rho,
        rho*nrg,
        rho * vel.x, rho * vel.y, rho * vel.z,
    };
    flux = {
        rho * vel_n,
        rho * vel_n * ent,
        rho * vel_n * vel.x + p * n.x,
        rho * vel_n * vel.y + p * n.y,
        rho * vel_n * vel.z + p * n.z,
    };
}

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

typedef class MhdHydroVars MhdFluidVarsIdealGas;
class MhdPhysicsIdealGas {
public:
    static constexpr int_t num_vars = 5;
    typedef MhdFluidVarsIdealGas MhdFluidStateT;
};

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //