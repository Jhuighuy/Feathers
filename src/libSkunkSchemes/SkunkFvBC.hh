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

#pragma once
#ifndef MHD_FV_BC_HH
#define MHD_FV_BC_HH

#include "SkunkBase.hh"
#include "libSkunkMesh/SkunkMesh.hh"
#include "SkunkHydro.hh"
#include <functional>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/**
 * Abstract finite-volume boundary condition.
 */
template<int_t num_vars_t>
class MhdFvBcT : public std::enable_shared_from_this<MhdFvBcT<num_vars_t>> {
public:
    /** Compute the ghost state values. */
    void get_ghost_state(TScalarField<num_vars_t>& u) const {
        get_ghost_state_(u);
    }
private:
    /** Compute the ghost state values. */
    virtual void get_ghost_state_(TScalarField<num_vars_t>& u) const = 0;
};  // class MhdFvBcT

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/**
 * @brief Abstract boundary condition.
 */
template<typename MhdPhysicsT>
class MhdFvBcPT :
    public std::enable_shared_from_this<MhdFvBcPT<MhdPhysicsT>> {
public:
    using MhdFluidStateT = typename MhdPhysicsT::MhdFluidStateT;
    static constexpr int_t num_vars = MhdPhysicsT::num_vars;

public:
    /** @brief Compute the ghost states. */
    void get_ghost_state(const mhd_vec3_t& n,
                         const mhd_vec3_t& r, const mhd_vec3_t& r_ghost,
                         const std::array<real_t, num_vars>& u,
                         std::array<real_t, num_vars>& u_ghost) const {
        MhdFluidStateT u_state(n, u.data()), u_ghost_state;
        get_ghost_state_(n, r, r_ghost, u_state, u_ghost_state);
        u_ghost_state.make_cons();
        u_ghost = u_ghost_state.cons;
    }

private:
    /** @brief Compute the ghost state. */
    virtual void get_ghost_state_(const mhd_vec3_t& n, const mhd_vec3_t& r, const mhd_vec3_t& r_ghost,
                                  const MhdFluidStateT& u,
                                  MhdFluidStateT& u_ghost) const = 0;
};  // MhdFvBcT

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/**
 * @brief Far field boundary condition.
 * Sets ghost state values to infinity state.
 */
template<typename MhdPhysicsT>
class MhdFvBcFarFieldT : public MhdFvBcPT<MhdPhysicsT> {
public:
    using typename MhdFvBcPT<MhdPhysicsT>::MhdFluidStateT;
    using MhdFvBcPT<MhdPhysicsT>::num_vars;

private:
    /** @brief Compute the ghost state. */
    void get_ghost_state_(const mhd_vec3_t& n, const mhd_vec3_t& r, const mhd_vec3_t& r_ghost,
                          const MhdFluidStateT& u,
                          MhdFluidStateT& u_ghost) const override;
};  // class MhdFvBcFarFieldT

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/**
 * @brief No-slip wall boundary condition.
 */
template<typename MhdPhysicsT>
class MhdFvBcNoSlipT : public MhdFvBcPT<MhdPhysicsT> {
public:
    using typename MhdFvBcPT<MhdPhysicsT>::MhdFluidStateT;
    using MhdFvBcPT<MhdPhysicsT>::num_vars;
    std::function<mhd_vec3_t(mhd_vec3_t)> vfunc;

    MhdFvBcNoSlipT(std::function<mhd_vec3_t(mhd_vec3_t)> vfunc = nullptr) : vfunc(std::move(vfunc)) {}

private:
    /** @brief Compute the ghost state. */
    void get_ghost_state_(const mhd_vec3_t& n, const mhd_vec3_t& r, const mhd_vec3_t& r_ghost,
                          const MhdFluidStateT& u,
                          MhdFluidStateT& u_ghost) const override;
};  // class MhdFvBcNoSlipT

/**
 * @brief Slip wall boundary condition.
 */
template<typename MhdPhysicsT>
class MhdFvBcSlipT :
    public MhdFvBcPT<MhdPhysicsT> {
public:
    using typename MhdFvBcPT<MhdPhysicsT>::MhdFluidStateT;
    using MhdFvBcPT<MhdPhysicsT>::num_vars;

private:
    /** @brief Compute the ghost state. */
    void get_ghost_state_(const mhd_vec3_t& n, const mhd_vec3_t& r, const mhd_vec3_t& r_ghost,
                          const MhdFluidStateT& u,
                          MhdFluidStateT& u_ghost) const override;
};  // class MhdFluidBcSlipTMhdFluidBcSlipT

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef MHD_FV_BC_HH
