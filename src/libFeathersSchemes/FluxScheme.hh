/*
 *  ______  ______   ______   ______  __  __   ______   ______   ______
 * /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\
 * \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \
 *  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\
 *   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
 *
 * Copyright (c) 2021 Oleg Butakov
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
#ifndef FLUX_SCHEME_HH_
#define FLUX_SCHEME_HH_

#include "SkunkBase.hh"
#include "SkunkHydro.hh"
//#include "SkunkFluidPhysics.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Abstract numerical flux.
 */
template<int_t num_vars_t>
class iFluxScheme : public tObject<iFluxScheme<num_vars_t>> {
public:
    /** Compute the numerical flux. */
    virtual void get_numerical_flux(const vec3_t& n,
                                    const std::array<real_t, num_vars_t>& ur,
                                    const std::array<real_t, num_vars_t>& ul,
                                    std::array<real_t, num_vars_t>& flux) const = 0;
}; // class iFluxScheme

/**
 * Abstract physics-based numerical flux.
 */
template<typename tPhysics>
class iPhysFluxScheme : public iFluxScheme<tPhysics::num_vars> {
public:
    static constexpr int_t num_vars = tPhysics::num_vars;
    using tFluidState = typename tPhysics::MhdFluidStateT;

public:
    /** Compute the numerical flux. */
    /** @{ */
    void get_numerical_flux(const vec3_t& n,
                            const std::array<real_t, num_vars>& ur,
                            const std::array<real_t, num_vars>& ul,
                            std::array<real_t, num_vars>& f) const final {
        get_numerical_flux(n,
                           tFluidState(n, ur.data()),
                           tFluidState(n, ul.data()), f);
    }
    virtual void get_numerical_flux(const vec3_t& n,
                                    const tFluidState& ur,
                                    const tFluidState& ul,
                                    std::array<real_t, num_vars>& f) const = 0;
    /** @} */
}; // class iPhysFluxScheme

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * @brief Local Lax-Friedrichs (Rusanov) numerical flux.
 *
 * Use this numerical flux if all other fails. 
 * It should always work.
 */
template<typename tPhysics>
class tLaxFriedrichsFluxScheme final : public iPhysFluxScheme<tPhysics> {
public:
    using iPhysFluxScheme<tPhysics>::num_vars;
    using typename iPhysFluxScheme<tPhysics>::tFluidState;

    /** Compute the numerical flux. */
    void get_numerical_flux(const vec3_t& n,
                            const tFluidState& ur,
                            const tFluidState& ul,
                            std::array<real_t, num_vars>& f) const final;
}; // class tLaxFriedrichsFluxScheme

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * @brief Harten-Lax-van Leer-Einfeldt numerical flux.
 *
 * Use this numerical flux if HLLC fails. 
 * It should (almost) always work.
 */
template<typename tPhysics>
class tHllFluxScheme : public iPhysFluxScheme<tPhysics> {
public:
    using iPhysFluxScheme<tPhysics>::num_vars;
    using typename iPhysFluxScheme<tPhysics>::tFluidState;

    /** Compute the signal speed. */
    void get_signal_speed(const tFluidState& ur,
                          const tFluidState& ul,
                          real_t& sr, real_t& sl) const;

    /** Compute the numerical flux. */
    void get_numerical_flux(const vec3_t& n,
                            const tFluidState& ur,
                            const tFluidState& ul,
                            std::array<real_t, num_vars>& f) const final;
}; // class tHllFluxScheme

/**
 * @brief Harten-Lax-van Leer-Contact numerical flux.
 *
 * Optimal choice for both gas and plasma physics.
 * In plasma physics case may be a bit more dissipative, but more consistent than HLLD/Roe.
 */
template<typename tPhysics>
class tHllcFluxScheme : public iPhysFluxScheme<tPhysics> {
public:
    using iPhysFluxScheme<tPhysics>::num_vars;
    using typename iPhysFluxScheme<tPhysics>::tFluidState;

    /** Compute the signal speed. */
    void get_signal_speed(const tFluidState& ur,
                          const tFluidState& ul,
                          real_t& sr, real_t& sl) const;

    /** Compute the numerical flux. */
    void get_numerical_flux(const vec3_t& n,
                            const tFluidState& ur,
                            const tFluidState& ul,
                            std::array<real_t, num_vars>& f) const override;
}; // class tHllcFluxScheme

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * @brief Roe numerical flux.
 *
 * Another optimal choice for gas physics.
 * For plasma physics is significantly slower that the HLLC/HLLD fluxes,
 * but sometimes may produce great results.
 */
template<typename tPhysics>
class tRoeFluxScheme : public iPhysFluxScheme<tPhysics> {
public:
    using iPhysFluxScheme<tPhysics>::num_vars;
    using typename iPhysFluxScheme<tPhysics>::tFluidState;

    /** Compute the numerical flux. */
    void get_numerical_flux_(const vec3_t& n,
                             const tFluidState& ur,
                             const tFluidState& ul,
                             std::array<real_t, num_vars>& f) const override;
}; // class tRoeFluxScheme

} // namespace feathers

#endif  // ifndef FLUX_SCHEME_HH_
