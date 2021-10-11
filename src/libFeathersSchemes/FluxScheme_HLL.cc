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

#include "FluxScheme.hh"

namespace {
const double gamma = Gamma;
const double gamma_2 = (gamma + 1.0)/(2.0*gamma);
}

namespace feathers {

/**
 * Calculate the Harten-Lax-van Leer-Einfeldt numerical flux.
 * @verbatim
 * [1] Eleuterio F. Toro,
 *     "Riemann Solvers and Numerical Methods
 *      for Fluid Dynamics" (Third Edition, 2009).
 * @endverbatim
 */
void tHllFluxScheme<tGasPhysics>::get_numerical_flux(const vec3_t& n,
                                                     const tFluidState& ur,
                                                     const tFluidState& ul,
                                                     real_t* f) const {
    /* Calculate Roe average sound speed.
     * [1] Eq. (10.53-10.54). */
    const real_t rr = std::sqrt(ur.rho);
    const real_t rl = std::sqrt(ul.rho);
    const real_t rs = 0.5*rr*rl/std::pow(rr + rl, 2);
    const real_t cs = std::sqrt(
        (rr*std::pow(ur.c_snd, 2) + rl*std::pow(ul.c_snd, 2))/(rr + rl) +
        rs*std::pow(ur.vel_n - ul.vel_n, 2));

    /* Calculate signal speeds.
     * [1], Eq. (10.52). */
    const real_t sr = 0.5*(ur.vel_n + ul.vel_n) + cs;
    const real_t sl = 0.5*(ur.vel_n + ul.vel_n) - cs;

    /* Supersonic cases.
     * [1], Eq. (10.20-10.21). */
    if (sr <= 0.0) {
        for (uint_t i = 0; i < num_vars; ++i) f[i] = ur.flux[i];
        return;
    }
    if (sl >= 0.0) {
        for (uint_t i = 0; i < num_vars; ++i) f[i] = ul.flux[i];
        return;
    }

    /* Subsonic case.
     * [1], Eq. (10.20-10.21). */
    if (sl <= 0.0 && 0.0 <= sr) {
        const real_t is = 1.0/(sr - sl);
        for (uint_t i = 0; i < num_vars; ++i) {
            f[i] = is*(sr*ul.flux[i] - sl*ur.flux[i] + sr*sl*(ur.cons[i] - ul.cons[i]));
        }
        return;
    }

    FEATHERS_ENSURE(!"Broken signal velocities.");
} // tHllFluxScheme::get_numerical_flux

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Calculate the Harten-Lax-van Leer-Contact numerical flux.
 * @verbatim
 * [1] Eleuterio F. Toro,
 *     "Riemann Solvers and Numerical Methods
 *      for Fluid Dynamics" (Third Edition, 2009).
 * @endverbatim
 */
void tHllcFluxScheme<tGasPhysics>::get_numerical_flux(const vec3_t& n,
                                                      const tFluidState& ur,
                                                      const tFluidState& ul,
                                                      real_t* f) const {
    /* Calculate average variables.
     * [1], Eq. (10.61-10.62). */
    const real_t rho = 0.5*(ur.rho + ul.rho);
    const real_t c_snd = 0.5*(ur.c_snd + ul.c_snd);
    const real_t p = std::max(0.0, 0.5*(ur.p + ul.p - rho*c_snd*(ur.vel_n - ul.vel_n)));

    /* Calculate sound speed coefficients.
     * [1], Eq. (10.60). */
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

    /* Calculate signal speeds.
     * [1], Eq. (10.59). */
    const real_t sr = ur.vel_n + ur.c_snd*gp;
    const real_t sl = ul.vel_n - ul.c_snd*gm;

    /* Supersonic cases.
     * [1], Eq. (10.20-10.21). */
    if (sr <= 0.0) {
        for (uint_t i = 0; i < num_vars; ++i) f[i] = ur.flux[i];
        return;
    }
    if (sl >= 0.0) {
        for (uint_t i = 0; i < num_vars; ++i) f[i] = ul.flux[i];
        return;
    }

    /* Subsonic cases.
     * [1], Eq. (10.37-10.39). */
    const real_t ss =
        ((ur.rho*ur.vel_n*(sr - ur.vel_n) - ur.p) -
         (ul.rho*ul.vel_n*(sl - ul.vel_n) - ul.p)) /
        (ur.rho*(sr - ur.vel_n) - ul.rho*(sl - ul.vel_n));
    if (ss <= 0.0 && 0.0 <= sr) {
        tFluidState us;
        const real_t is = 1.0/(sr - ss);
        us.rho = ur.rho*(sr - ur.vel_n)*is;
        us.nrg = ur.nrg + (ss - ur.vel_n)*(ss + ur.p/ur.rho*is);
        us.vel = ur.vel + (ss - ur.vel_n)*n;
        us.make_cons();
        for (uint_t i = 0; i < num_vars; ++i) {
            f[i] = ur.flux[i] + sr*(us.cons[i] - ur.cons[i]);
        }
        return;
    }
    if (sl <= 0.0 && 0.0 <= ss)  {
        tFluidState us;
        const real_t is = 1.0/(sl - ss);
        us.rho = ul.rho*(sl - ul.vel_n)*is;
        us.nrg = ul.nrg + (ss - ul.vel_n)*(ss + ul.p/ul.rho*is);
        us.vel = ul.vel + (ss - ul.vel_n)*n;
        us.make_cons();
        for (uint_t i = 0; i < num_vars; ++i) {
            f[i] = ul.flux[i] + sl*(us.cons[i] - ul.cons[i]);
        }
        return;
    }

    FEATHERS_ENSURE(!"Broken signal velocities.");
} // tHllcFluxScheme<tGasPhysics>::get_numerical_flux

} // namespace feathers
