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

#include "FluxScheme.hh"

namespace {
const double gamma = Gamma;
const double gamma_2 = (gamma + 1.0)/(2.0*gamma);
}

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * @brief Calculate signal speed estimates.
 * @verbatim
 * [1] Eleuterio F. Toro,
 *     "Riemann Solvers and Numerical Methods
 *      for Fluid Dynamics" (Third Edition, 2009).
 * @endverbatim
 */
template<>
void tHLLFluxScheme<MhdPhysicsIdealGas>::get_signal_speed(const tFluidState& ur,
                                                          const tFluidState& ul,
                                                          real_t& sr, real_t& sl) const {
    /*
     * Calculate Roe average sound speed.
     * [1] Eq. (10.53-10.54).
     */
    const real_t rr = std::sqrt(ur.rho);
    const real_t rl = std::sqrt(ul.rho);
    const real_t rs = 0.5 * rr * rl / std::pow(rr + rl, 2);
    const real_t cs = std::sqrt((rr * ur.c2snd + rl * ul.c2snd) / (rr + rl) + rs*std::pow(ur.Vn - ul.Vn, 2));
    /*
     * Calculate signal speeds.
     * [1], Eq. (10.52).
     */
    sr = 0.5*(ur.Vn + ul.Vn) + cs;
    sl = 0.5*(ur.Vn + ul.Vn) - cs;
}   // tHLLFluxScheme<MhdPhysicsIdealGas>::get_signal_speed_

/**
 * @brief Calculate the Harten-Lax-van Leer-Einfeldt numerical flux.
 * @verbatim
 * [1] Eleuterio F. Toro,
 *     "Riemann Solvers and Numerical Methods
 *      for Fluid Dynamics" (Third Edition, 2009).
 * @endverbatim
 */
template<typename TPhysics>
void tHLLFluxScheme<TPhysics>::get_numerical_flux(const vec3_t& n,
                                                  const tFluidState& ur,
                                                  const tFluidState& ul,
                                                  std::array<real_t, num_vars>& f) const {
    /*
     * Supersonic cases.
     * [1], Eq. (10.20-10.21).
     */
    real_t sr, sl;
    get_signal_speed(ur, ul, sr, sl);
    if (sr <= 0.0) {
        f = ur.flux;
        return;
    }
    if (sl >= 0.0) {
        f = ul.flux;
        return;
    }
    /*
     * Subsonic case.
     * [1], Eq. (10.20-10.21).
     */
    if (sl <= 0.0 && 0.0 <= sr) {
        const real_t is = 1.0 / (sr - sl);
        for (int_t i = 0; i < num_vars; ++i) {
            f[i] = is*(sr*ul.flux[i] - sl*ur.flux[i] + sr*sl*(ur.cons[i] - ul.cons[i]));
        }
        return;
    }
    FEATHERS_ENSURE(!"Broken signal velocities.");
}   // tHLLFluxScheme::get_numerical_flux

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * @brief Calculate pressure-based signal speed estimates.
 * @verbatim
 * [1] Eleuterio F. Toro,
 *     "Riemann Solvers and Numerical Methods
 *      for Fluid Dynamics" (Third Edition, 2009).
 * @endverbatim
 */
template<>
void tHLLCFluxScheme<MhdPhysicsIdealGas>::get_signal_speed(const tFluidState& ur,
                                                           const tFluidState& ul,
                                                           real_t& sr, real_t& sl) const {
    /*
     * Calculate average variables.
     * [1], Eq. (10.61-10.62).
     */
    const real_t rho = 0.5 * (ur.rho + ur.rho);
    const real_t c_snd = 0.5 * (ur.c_snd + ur.c_snd);
    const real_t p = std::max(0.5 * (ur.p + ul.p - rho * c_snd * (ur.Vn - ul.Vn)), 0.0);
    /*
     * Calculate sound speed coefficients.
     * [1], Eq. (10.60).
     */
    real_t gp;
    if (p > ur.p) {
        gp = std::sqrt(1.0 + /*m_phys->*/gamma_2*(p/ur.p - 1.0));
    } else {
        gp = 1.0;
    }
    real_t gm;
    if (p > ul.p) {
        gm = std::sqrt(1.0 + /*m_phys->*/gamma_2*(p/ul.p - 1.0));
    } else {
        gm = 1.0;
    }
    /*
     * Calculate signal speeds.
     * [1], Eq. (10.59).
     */
    sr = ur.Vn + ur.c_snd*gp;
    sl = ul.Vn - ul.c_snd*gm;
}   // tHLLCFluxScheme<MhdPhysicsIdealGas>::get_signal_speed_

#define HLLC_VARIATION 0

/**
 * @brief Calculate the Harten-Lax-van Leer-Contact numerical flux.
 * @verbatim
 * [1] Eleuterio F. Toro,
 *     "Riemann Solvers and Numerical Methods
 *      for Fluid Dynamics" (Third Edition, 2009).
 * @endverbatim
 */
template<>
void tHLLCFluxScheme<MhdPhysicsIdealGas>::get_numerical_flux(const vec3_t& n,
                                                             const tFluidState& ur,
                                                             const tFluidState& ul,
                                                             std::array<real_t, num_vars>& f) const {
    /*
     * Supersonic cases.
     * [1], Eq. (10.20-10.21).
     */
    real_t sr, sl;
    get_signal_speed(ur, ul, sr, sl);
    if (sr <= 0.0) {
        f = ur.flux;
        return;
    }
    if (sl >= 0.0) {
        f = ul.flux;
        return;
    }
    /*
     * Subsonic cases.
     * [1], Eq. (10.37).
     */
    const real_t ss =
        ((ur.rho*ur.Vn*(sr - ur.Vn) - ur.p) -
         (ul.rho*ul.Vn*(sl - ul.Vn) - ul.p))/(ur.rho*(sr - ur.Vn) - ul.rho*(sl - ul.Vn));
#if HLLC_VARIATION == 0
    /*
     * Original HLLC:
     * [1], Eq. (10.38), (10.39).
     */
    if (ss <= 0.0 && 0.0 <= sr) {
        tFluidState us;
        us.rho = ur.rho*(sr - ur.Vn)/(sr - ss);
        us.nrg = ur.nrg + (ss - ur.Vn)*(ss + ur.p/ur.rho/(sr - ss));
        us.V   = ur.V - n*(ur.Vn - ss);
        us.make_cons();
        for (int_t i = 0; i < num_vars; ++i) {
            f[i] = ur.flux[i] + sr*(us.cons[i] - ur.cons[i]);
        }
        return;
    }
    if (sl <= 0.0 && 0.0 <= ss)  {
        tFluidState us;
        us.rho = ul.rho*(sl - ul.Vn)/(sl - ss);
        us.nrg = ul.nrg + (ss - ul.Vn)*(ss + ul.p/ul.rho/(sl - ss));
        us.V   = ul.V - n*(ul.Vn - ss);
        us.make_cons();
        for (int_t i = 0; i < num_vars; ++i) {
            f[i] = ul.flux[i] + sl*(us.cons[i] - ul.cons[i]);
        }
        return;
    }
#elif HLLC_VARIATION == 1
    /*
     * Variation 1 of HLLC:
     * [1], Eq. (10.40), (10.41).
     */
    const std::array<real_t, 5> ds{ 0.0, ss, n.x, n.y, n.z };
    if (ss <= 0.0 && 0.0 <= sr) {
        const real_t is = 1.0/(sr - ss);
        for (int_t i = 0; i < num_vars; ++i) {
            f[i] = is*(ss*(sr*ur.cons[i] - ur.flux[i]) +
                       sr*(ur.p + ur.rho*(sr - ur.Vn)*(ss - ur.Vn))*ds[i]);
        }
        return;
    }
    if (sl <= 0.0 && 0.0 <= ss) {
        const real_t is = 1.0/(sl - ss);
        for (int_t i = 0; i < num_vars; ++i) {
            f[i] = is*(ss*(sl*ul.cons[i] - ul.flux[i]) +
                       sl*(ul.p + ul.rho*(sl - ul.Vn)*(ss - ul.Vn))*ds[i]);
        }
        return;
    }
#elif HLLC_VARIATION == 2
    /*
     * Variation 2 of HLLC:
     * [1], Eq. (10.40), (10.42), (10.44).
     */
    const std::array<real_t, 5> ds{ 0.0, ss, n.x, n.y, n.z };
    const real_t ps = 0.5*((ur.p + ur.rho*(sr - ur.Vn)*(ss - ur.Vn)) +
                               (ul.p + ul.rho*(sl - ul.Vn)*(ss - ul.Vn)));
    if (ss <= 0.0 && 0.0 <= sr) {
        const real_t is = 1.0/(sr - ss);
        for (int_t i = 0; i < num_vars; ++i) {
            f[i] = is*(ss*(sr*ur.cons[i] - ur.flux[i]) + sr*ps*ds[i]);
        }
        return;
    }
    if (sl <= 0.0 && 0.0 <= ss) {
        const real_t is = 1.0/(sl - ss);
        for (int_t i = 0; i < num_vars; ++i) {
            f[i] = is*(ss*(sl*ul.cons[i] - ul.flux[i]) + sl*ps*ds[i]);
        }
        return;
    }
#endif
    assert(!"Broken signal velocities.");
}   // tHLLCFluxScheme<MhdPhysicsIdealGas>::get_numerical_flux

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template class feathers::tHLLFluxScheme<MhdPhysicsIdealGas>;
template class feathers::tHLLCFluxScheme<MhdPhysicsIdealGas>;
