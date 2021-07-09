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
#ifndef FLUX_SCHEME_HH_
#define FLUX_SCHEME_HH_

#include "SkunkBase.hh"
#include "SkunkHydro.hh"
//#include "SkunkFluidPhysics.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/**
 * Abstract numerical flux.
 */
template<int_t num_vars_t>
class IFluxScheme : public TObject<IFluxScheme<num_vars_t>> {
public:
    /** Compute the numerical flux. */
    virtual void get_numerical_flux(const vec3_t& n,
                                    const std::array<real_t, num_vars_t>& ur,
                                    const std::array<real_t, num_vars_t>& ul,
                                    std::array<real_t, num_vars_t>& flux) const = 0;
};  // class IFluxScheme

/**
 * Abstract physics-based numerical flux.
 */
template<typename TPhysics>
class IPhysicalFluxScheme : public IFluxScheme<TPhysics::num_vars> {
public:
    static constexpr int_t num_vars = TPhysics::num_vars;
    using TFluidState = typename TPhysics::MhdFluidStateT;

public:
    /** Compute the numerical flux. */
    /** @{ */
    void get_numerical_flux(const vec3_t& n,
                            const std::array<real_t, num_vars>& ur,
                            const std::array<real_t, num_vars>& ul,
                            std::array<real_t, num_vars>& f) const final {
        get_numerical_flux(n,
                           TFluidState(n, ur.data()),
                           TFluidState(n, ul.data()), f);
    }
    virtual void get_numerical_flux(const vec3_t& n,
                                    const TFluidState& ur,
                                    const TFluidState& ul,
                                    std::array<real_t, num_vars>& f) const = 0;
    /** @} */
};  // class IPhysicalFluxScheme

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/**
 * @brief Local Lax-Friedrichs (Rusanov) numerical flux.
 *
 * Use this numerical flux if all other fails. 
 * It should always work.
 */
template<typename TPhysics>
class TLaxFriedrichsFluxScheme final : public IPhysicalFluxScheme<TPhysics> {
public:
    using IPhysicalFluxScheme<TPhysics>::num_vars;
    using typename IPhysicalFluxScheme<TPhysics>::TFluidState;

    /** Compute the numerical flux. */
    void get_numerical_flux(const vec3_t& n,
                            const TFluidState& ur,
                            const TFluidState& ul,
                            std::array<real_t, num_vars>& f) const final;
};  // class TLaxFriedrichsFluxScheme

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/**
 * @brief Harten-Lax-van Leer-Einfeldt numerical flux.
 *
 * Use this numerical flux if HLLC fails. 
 * It should (almost) always work.
 */
template<typename TPhysics>
class THLLFluxScheme : public IPhysicalFluxScheme<TPhysics> {
public:
    using IPhysicalFluxScheme<TPhysics>::num_vars;
    using typename IPhysicalFluxScheme<TPhysics>::TFluidState;

    /** Compute the signal speed. */
    void get_signal_speed(const TFluidState& ur,
                          const TFluidState& ul,
                          real_t& sr, real_t& sl) const;

    /** Compute the numerical flux. */
    void get_numerical_flux(const vec3_t& n,
                            const TFluidState& ur,
                            const TFluidState& ul,
                            std::array<real_t, num_vars>& f) const final;
};  // class THLLFluxScheme

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * @brief Harten-Lax-van Leer-Contact numerical flux.
 *
 * Optimal choice for both gas and plasma physics.
 * In plasma physics case may be a bit more dissipative, but more consistent than HLLD/Roe.
 */
template<typename TPhysics>
class THLLCFluxScheme : public IPhysicalFluxScheme<TPhysics> {
public:
    using IPhysicalFluxScheme<TPhysics>::num_vars;
    using typename IPhysicalFluxScheme<TPhysics>::TFluidState;

    /** Compute the signal speed. */
    void get_signal_speed(const TFluidState& ur,
                          const TFluidState& ul,
                          real_t& sr, real_t& sl) const;

    /** Compute the numerical flux. */
    void get_numerical_flux(const vec3_t& n,
                            const TFluidState& ur,
                            const TFluidState& ul,
                            std::array<real_t, num_vars>& f) const override;
};  // class THLLCFluxScheme

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/**
 * @brief Roe numerical flux.
 *
 * Another optimal choice for gas physics.
 * For plasma physics is significantly slower that the HLLC/HLLD fluxes,
 * but sometimes may produce great results.
 */
template<typename TPhysics>
class TRoeFluxScheme : public IPhysicalFluxScheme<TPhysics> {
public:
    using IPhysicalFluxScheme<TPhysics>::num_vars;
    using typename IPhysicalFluxScheme<TPhysics>::TFluidState;

    /** Compute the numerical flux. */
    void get_numerical_flux_(const vec3_t& n,
                             const TFluidState& ur,
                             const TFluidState& ul,
                             std::array<real_t, num_vars>& f) const override;
};  // class TRoeFluxScheme

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef FLUX_SCHEME_HH_
