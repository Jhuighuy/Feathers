// Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver.
// Copyright(C) Butakov Oleg and Co. 2019.

#pragma once

#include "SkunkBase.hh"
#include "libSkunkMesh/SkunkMesh.hh"
#include <cstring>

using mhd_vec3_t = skunk::vec3_t;

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

static const real_t Gamma = 1.4;
static const real_t Gamma1 = Gamma - 1.0;

class MhdHydroVars {
public:
    real_t rho;
    union {
        struct { real_t nrg, eps, kin; };
        struct { real_t e, e_int, e_kin; };
    };
    real_t ent, p;
    union {
        mhd_vec3_t V;
        mhd_vec3_t vel;
    };
    union {
        struct { real_t Vn, V2; };
        struct { real_t vel_n, vel_sqr; };
    };
    real_t c_snd, c2snd;
    std::array<real_t, 5> prim;
    std::array<real_t, 5> cons;
    std::array<real_t, 5> flux;

public:
    explicit MhdHydroVars() {
        std::memset((void*)this, 0, sizeof(*this));
    }
    explicit MhdHydroVars(const mhd_vec3_t& n,
                          const real_t* q_cons,
                          const real_t* q_prim = nullptr);
public:
    void make_cons() {
        //real_t pr[]{rho, p, V.x, V.y, V.z};
        //cons = MhdHydroVars({}, nullptr, pr).cons;
#if 1
        cons = {
                rho,
                rho*nrg,
                rho*V.x, rho*V.y, rho*V.z,
        };
#endif
    }

    void make_roe_average(const mhd_vec3_t& n,
                          const MhdHydroVars& qr, const MhdHydroVars& ql);

    skunk::vec5_t make_eigenvalues(const mhd_vec3_t& n) const;
    skunk::vec5_t make_roe_eigenvalues(const mhd_vec3_t& n) const;
    skunk::vec5_t make_roe_eigenvalues(const mhd_vec3_t& n,
                                       const MhdHydroVars& qr, const MhdHydroVars& ql) const;

    std::pair<skunk::mat5_t, skunk::mat5_t> make_eigenvectors(const mhd_vec3_t& n) const;
};

inline
MhdHydroVars::MhdHydroVars(const mhd_vec3_t& n,
                           const real_t* q_cons, const real_t* q_prim) : MhdHydroVars() {
    if (q_cons != nullptr) {
        rho = q_cons[0];
        if (rho == 0.0) return;
        nrg = q_cons[1]/rho;
        V.x = q_cons[2]/rho;
        V.y = q_cons[3]/rho;
        V.z = q_cons[4]/rho;
        Vn  = V.inner_prod(n);
        V2  = V.inner_prod(V);
        kin = 0.5*V2;
        eps = nrg - kin;
        p   = Gamma1*rho*eps;
        ent = nrg + p/rho;
    } else if (q_prim != nullptr) {
        rho = q_prim[0];
        p   = q_prim[1];
        V.x = q_prim[2];
        V.y = q_prim[3];
        V.z = q_prim[4];
        Vn  = V.inner_prod(n);
        V2  = V.inner_prod(V);
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
        rho, p, V.x, V.y, V.z,
    };
    cons = {
        rho,
        rho*nrg,
        rho*V.x, rho*V.y, rho*V.z,
    };
    flux = {
        rho*Vn,
        rho*Vn*ent,
        rho*Vn*V.x + p*n.x,
        rho*Vn*V.y + p*n.y,
        rho*Vn*V.z + p*n.z,
    };
}

inline
void MhdHydroVars::make_roe_average(const mhd_vec3_t& n,
                                    const MhdHydroVars& qr, const MhdHydroVars& ql) {
    const real_t rr = std::sqrt(qr.rho),
                     rl = std::sqrt(ql.rho);
    rho = rr*rl;
    ent = (rr*qr.ent + rl*ql.ent)/(rr + rl);
    vel = (rr*qr.vel + rl*ql.vel)/(rr + rl);
    vel_n = vel.inner_prod(n);
    e_kin = 0.5*vel.inner_prod(vel);
    //c2snd = (rr*qr.c2snd + rl*ql.c2snd)/(rr + rl)
    //      + 0.5*Gamma1*rr*rl*std::pow((qr.vel_n - ql.vel_n)/(rr + rl), 2);
    c2snd = Gamma1*(ent - e_kin);
    c_snd = std::sqrt(c2snd);
}   // MhdHydroVars::make_roe_average

inline
skunk::vec5_t MhdHydroVars::make_eigenvalues(const mhd_vec3_t& n) const {
    using namespace skunk;
    const vec5_t eig_vals {
        vel_n - c_snd,
        vel_n,
        vel_n + c_snd,
        vel_n,
        vel_n,
    };
    return eig_vals;
}   // MhdHydroVars::make_eigenvalues
inline
skunk::vec5_t MhdHydroVars::make_roe_eigenvalues(const mhd_vec3_t& n) const {
    using namespace skunk;
    const real_t abs_vel_n = std::abs(vel_n);
    const vec5_t eig_vals {
        std::abs(vel_n - c_snd),
        abs_vel_n,
        std::abs(vel_n + c_snd),
        abs_vel_n,
        abs_vel_n,
    };
    return eig_vals;
}   // MhdHydroVars::make_eigenvalues
inline
skunk::vec5_t MhdHydroVars::make_roe_eigenvalues(const mhd_vec3_t& n,
                                                 const MhdHydroVars& qr, const MhdHydroVars& ql) const {
    using namespace skunk;
    const real_t abs_vel_n = std::abs(vel_n);
    const vec5_t eig_vals {
        std::min(std::abs(vel_n - c_snd), std::abs(ql.vel_n - ql.c_snd)),
        abs_vel_n,
        std::max(std::abs(vel_n + c_snd), std::abs(qr.vel_n + qr.c_snd)),
        abs_vel_n,
        abs_vel_n,
    };
    return eig_vals;
}   // MhdHydroVars::make_eigenvalues

inline
std::pair<skunk::mat5_t, skunk::mat5_t> MhdHydroVars::make_eigenvectors(const mhd_vec3_t& n) const {
    using namespace skunk;
    /* Compute the normal-independent parts of the Eigenvectors. */
    const vec3_t vel1 = Gamma1*vel;
    const real_t e_kin1 = Gamma1*e_kin;
    const mat5_t eig_vec_r_1{
        1.0, 1.0, 1.0, 0.0, 0.0,
        vel.x - c_snd*n.x, vel.x, vel.x + c_snd*n.x, 0.0, 0.0,
        vel.y - c_snd*n.y, vel.y, vel.y + c_snd*n.y, 0.0, 0.0,
        vel.z - c_snd*n.z, vel.z, vel.z + c_snd*n.z, 0.0, 0.0,
        ent - c_snd*vel_n, e_kin, ent + c_snd*vel_n, 0.0, 0.0,
    };
    const mat5_t eig_vec_l_1 = 0.5/c2snd*mat5_t{
        e_kin1 + c_snd*vel_n, -vel1.x - c_snd*n.x, -vel1.y - c_snd*n.y, -vel1.z - c_snd*n.z, Gamma1,
        2.0*(-e_kin1 + c2snd), 2.0*vel1.x, 2.0*vel1.y, 2.0*vel1.z, -2.0*Gamma1,
        e_kin1 - c_snd*vel_n, -vel1.x + c_snd*n.x, -vel1.y + c_snd*n.y, -vel1.z + c_snd*n.z, Gamma1,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
    };
    /* Compute the normal-dependent parts of the Eigenvectors. */
    mat5_t eig_vec_l_2, eig_vec_r_2;
    if (std::abs(n.x) > std::abs(n.y) && std::abs(n.x) > std::abs(n.z)) {
        SKUNK_ASSERT(n.x != 0.0);
        eig_vec_r_2 = {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, +n.y, -n.z,
            0.0, 0.0, 0.0, -n.x,  0.0,
            0.0, 0.0, 0.0,  0.0, +n.x,
            0.0, 0.0, 0.0, vel.x*n.y - vel.y*n.x, vel.z*n.x - vel.x*n.z,
        };
        eig_vec_l_2 = {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            +(vel.y - vel_n*n.y)/n.x, +n.y, -(1.0 - n.y*n.y)/n.x, +n.y*n.z/n.x, 0.0,
            -(vel.z - vel_n*n.z)/n.x, -n.z, -n.y*n.z/n.x, +(1.0 - n.z*n.z)/n.x, 0.0,
        };
    } else if (std::abs(n.y) > std::abs(n.x) && std::abs(n.y) > std::abs(n.z)) {
        SKUNK_ASSERT(n.y != 0.0);
        eig_vec_r_2 = {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, +n.y,  0.0,
            0.0, 0.0, 0.0, -n.x, +n.z,
            0.0, 0.0, 0.0,  0.0, -n.y,
            0.0, 0.0, 0.0, vel.x*n.y - vel.y*n.x, vel.y*n.z - vel.z*n.y,
        };
        eig_vec_l_2 = {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            -(vel.x - vel_n*n.x)/n.y, +(1.0 - n.x*n.x)/n.y, -n.x, +n.x*n.z/n.y, 0.0,
            +(vel.z - vel_n*n.z)/n.y, -n.x*n.z/n.y, +n.z, -(1.0 - n.z*n.z)/n.y, 0.0,
        };
    } else if (std::abs(n.z) > std::abs(n.x) && std::abs(n.z) > std::abs(n.y)) {
        SKUNK_ASSERT(n.z != 0.0);
        eig_vec_r_2 = {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -n.z,  0.0,
            0.0, 0.0, 0.0,  0.0, +n.z,
            0.0, 0.0, 0.0, +n.x, -n.y,
            0.0, 0.0, 0.0, vel.z*n.x - vel.x*n.z, vel.y*n.z - vel.z*n.y,
        };
        eig_vec_l_2 = {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            +(vel.x - vel_n*n.x)/n.z, -(1.0 - n.x*n.x)/n.z, +n.x*n.y/n.z, +n.x, 0.0,
            -(vel.y - vel_n*n.y)/n.z, -n.x*n.y/n.z, +(1.0 - n.y*n.y)/n.z, -n.y, 0.0,
        };
    }
    /* Gather the final results:
     * multiply by the permutation matrix. */
    static const mat5_t permute_direct {
        1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0,
        0.0, 1.0, 0.0, 0.0, 0.0,
    };
    static const mat5_t permute_reverse {
        1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0,
        0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0,
    };
    const mat5_t ev_r = permute_reverse.prod(eig_vec_r_1 + eig_vec_r_2);
    const mat5_t ev_l = (eig_vec_l_1 + eig_vec_l_2).prod(permute_direct);
    return std::make_pair(ev_r, ev_l);
}   // MhdHydroVars::make_eigenvectors

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