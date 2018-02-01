// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::PengRobinson
 */
#ifndef OPM_PENG_ROBINSON_HPP
#define OPM_PENG_ROBINSON_HPP

#include <opm/material/fluidstates/TemperatureOverlayFluidState.hpp>
#include <opm/material/IdealGas.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>

#include <opm/material/common/Unused.hpp>
#include <opm/material/common/PolynomialUtils.hpp>

#include <csignal>

namespace Opm {

/*!
 * \brief Implements the Peng-Robinson equation of state for liquids
 *        and gases.
 *
 * See:
 *
 * D.-Y. Peng, D.B. Robinson: A new two-constant equation of state,
 * Industrial & Engineering Chemistry Fundamentals, 1976, 15 (1),
 * pp. 59–64
 *
 * R. Reid, et al.: The Properties of Gases and Liquids,
 * 4th edition, McGraw-Hill, 1987, pp. 42-44, 82
 */
template <class Scalar>
class PengRobinson
{
    //! The ideal gas constant [Pa * m^3/mol/K]
    static const Scalar R;

    PengRobinson()
    { }

public:
    static void  init(Scalar /*aMin*/, Scalar /*aMax*/, unsigned /*na*/,
                      Scalar /*bMin*/, Scalar /*bMax*/, unsigned /*nb*/)
    {
/*
        // resize the tabulation for the critical points
        criticalTemperature_.resize(aMin, aMax, na, bMin, bMax, nb);
        criticalPressure_.resize(aMin, aMax, na, bMin, bMax, nb);
        criticalMolarVolume_.resize(aMin, aMax, na, bMin, bMax, nb);

        Scalar VmCrit, pCrit, TCrit;
        for (unsigned i = 0; i < na; ++i) {
            Scalar a = criticalTemperature_.iToX(i);
            assert(std::abs(criticalTemperature_.xToI(criticalTemperature_.iToX(i)) - i) < 1e-10);

            for (unsigned j = 0; j < nb; ++j) {
                Scalar b = criticalTemperature_.jToY(j);
                assert(std::abs(criticalTemperature_.yToJ(criticalTemperature_.jToY(j)) - j) < 1e-10);

                findCriticalPoint_(TCrit, pCrit, VmCrit,  a, b);

                criticalTemperature_.setSamplePoint(i, j, TCrit);
                criticalPressure_.setSamplePoint(i, j, pCrit);
                criticalMolarVolume_.setSamplePoint(i, j, VmCrit);
            }
        }
        */
    }

    /*!
     * \brief Predicts the vapor pressure for the temperature given in
     *        setTP().
     *
     * Initially, the vapor pressure is roughly estimated by using the
     * Ambrose-Walton method, then the Newton method is used to make
     * difference between the gas and liquid phase fugacity zero.
     */
    template <class Evaluation, class Params>
    static Evaluation computeVaporPressure(const Params& params, const Evaluation& T)
    {
        typedef typename Params::Component Component;
        if (T >= Component::criticalTemperature())
            return Component::criticalPressure();

        // initial guess of the vapor pressure
        Evaluation Vm[3];
        const Scalar eps = Component::criticalPressure()*1e-10;

        // use the Ambrose-Walton method to get an initial guess of
        // the vapor pressure
        Evaluation pVap = ambroseWalton_(params, T);

        // Newton-Raphson method
        for (unsigned i = 0; i < 5; ++i) {
            // calculate the molar densities
            int numSol OPM_OPTIM_UNUSED = molarVolumes(Vm, params, T, pVap);
            assert(numSol == 3);

            const Evaluation& f = fugacityDifference_(params, T, pVap, Vm[0], Vm[2]);
            Evaluation df_dp =
                fugacityDifference_(params, T, pVap  + eps, Vm[0], Vm[2])
                -
                fugacityDifference_(params, T, pVap - eps, Vm[0], Vm[2]);
            df_dp /= 2*eps;

            const Evaluation& delta = f/df_dp;
            pVap = pVap - delta;

            if (std::abs(Opm::scalarValue(delta/pVap)) < 1e-10)
                break;
        }

        return pVap;
    }

    /*!
     * \brief Computes molar volumes where the Peng-Robinson EOS is
     *        true.
     */
    template <class FluidState, class Params>
    static
    typename FluidState::Scalar
    computeMolarVolume(const FluidState& fs,
                       Params& params,
                       unsigned phaseIdx,
                       bool isGasPhase)
    {
        Valgrind::CheckDefined(fs.temperature(phaseIdx));
        Valgrind::CheckDefined(fs.pressure(phaseIdx));

        typedef typename FluidState::Scalar Evaluation;

        Evaluation Vm = 0;
        Valgrind::SetUndefined(Vm);

        const Evaluation& T = fs.temperature(phaseIdx);
        const Evaluation& p = fs.pressure(phaseIdx);

        const Evaluation& a = params.a(phaseIdx); // "attractive factor"
        const Evaluation& b = params.b(phaseIdx); // "co-volume"

        if (!std::isfinite(Opm::scalarValue(a))
            || std::abs(Opm::scalarValue(a)) < 1e-30)
            return std::numeric_limits<Scalar>::quiet_NaN();
        if (!std::isfinite(Opm::scalarValue(b)) || b <= 0)
            return std::numeric_limits<Scalar>::quiet_NaN();

        const Evaluation& RT= R*T;
        const Evaluation& Astar = a*p/(RT*RT);
        const Evaluation& Bstar = b*p/RT;

        const Evaluation& a1 = 1.0;
        const Evaluation& a2 = - (1 - Bstar);
        const Evaluation& a3 = Astar - Bstar*(3*Bstar + 2);
        const Evaluation& a4 = Bstar*(- Astar + Bstar*(1 + Bstar));

        // ignore the first two results if the smallest
        // compressibility factor is <= 0.0. (this means that if we
        // would get negative molar volumes for the liquid phase, we
        // consider the liquid phase non-existant.)
        Evaluation Z[3] = {0.0,0.0,0.0};
        Valgrind::CheckDefined(a1);
        Valgrind::CheckDefined(a2);
        Valgrind::CheckDefined(a3);
        Valgrind::CheckDefined(a4);
        int numSol = invertCubicPolynomial(Z, a1, a2, a3, a4);
        if (numSol == 3) {
            // the EOS has three intersections with the pressure,
            // i.e. the molar volume of gas is the largest one and the
            // molar volume of liquid is the smallest one
            if (isGasPhase)
                Vm = Z[2]*RT/p;
            else
                Vm = Z[0]*RT/p;
        }
        else if (numSol == 1) {
            // the EOS only has one intersection with the pressure,
            // for the other phase, we take the extremum of the EOS
            // with the largest distance from the intersection.
            Evaluation VmCubic = Z[0]*RT/p;
            Vm = VmCubic;

            // find the extrema (if they are present)
            Evaluation Vmin, Vmax, pmin, pmax;
            if (findExtrema_(Vmin, Vmax,
                             pmin, pmax,
                             a, b, T))
            {
                if (isGasPhase)
                    Vm = std::max(Vmax, VmCubic);
                else {
                    if (Vmin > 0)
                        Vm = std::min(Vmin, VmCubic);
                    else
                        Vm = VmCubic;
                }
            }
            else {
                // the EOS does not exhibit any physically meaningful
                // extrema, and the fluid is critical...
                Vm = VmCubic;
                handleCriticalFluid_(Vm, fs, params, phaseIdx, isGasPhase);
            }
        }

        Valgrind::CheckDefined(Vm);
        assert(std::isfinite(Opm::scalarValue(Vm)));
        assert(Vm > 0);
        return Vm;
    }

    /*!
     * \brief Returns the fugacity coefficient for a given pressure
     *        and molar volume.
     *
     * This is the same value as computeFugacity() because the mole
     * fraction of a component in a pure fluid is obviously always
     * 100%.
     *
     * \param params Parameters
     */
    template <class Evaluation, class Params>
    static Evaluation computeFugacityCoeffient(const Params& params)
    {
        const Evaluation& T = params.temperature();
        const Evaluation& p = params.pressure();
        const Evaluation& Vm = params.molarVolume();

        const Evaluation& RT = R*T;
        const Evaluation& Z = p*Vm/RT;
        const Evaluation& Bstar = p*params.b() / RT;

        const Evaluation& tmp =
            (Vm + params.b()*(1 + std::sqrt(2))) /
            (Vm + params.b()*(1 - std::sqrt(2)));
        const Evaluation& expo = - params.a()/(RT * 2 * params.b() * std::sqrt(2));
        const Evaluation& fugCoeff =
            Opm::exp(Z - 1) / (Z - Bstar) *
            Opm::pow(tmp, expo);

        return fugCoeff;
    }

    /*!
     * \brief Returns the fugacity coefficient for a given pressure
     *        and molar volume.
     *
     * This is the fugacity coefficient times the pressure. The mole
     * fraction of a component in a pure fluid is obviously always
     * 100%, so it is not required.
     *
     * \param params Parameters
     */
    template <class Evaluation, class Params>
    static Evaluation computeFugacity(const Params& params)
    { return params.pressure()*computeFugacityCoeff(params); }

protected:
    template <class FluidState, class Params, class Evaluation = typename FluidState::Scalar>
    static void handleCriticalFluid_(Evaluation& Vm,
                                     const FluidState& /*fs*/,
                                     const Params& params,
                                     unsigned phaseIdx,
                                     bool isGasPhase)
    {
        Evaluation Tcrit, pcrit, Vcrit;
        findCriticalPoint_(Tcrit,
                           pcrit,
                           Vcrit,
                           params.a(phaseIdx),
                           params.b(phaseIdx));


        //Evaluation Vcrit = criticalMolarVolume_.eval(params.a(phaseIdx), params.b(phaseIdx));

        if (isGasPhase)
            Vm = Opm::max(Vm, Vcrit);
        else
            Vm = Opm::min(Vm, Vcrit);
    }

    template <class Evaluation>
    static void findCriticalPoint_(Evaluation& Tcrit,
                                   Evaluation& pcrit,
                                   Evaluation& Vcrit,
                                   const Evaluation& a,
                                   const Evaluation& b)
    {
        Evaluation minVm(0);
        Evaluation maxVm(1e30);

        Evaluation minP(0);
        Evaluation maxP(1e30);

        // first, we need to find an isotherm where the EOS exhibits
        // a maximum and a minimum
        Evaluation Tmin = 250; // [K]
        for (unsigned i = 0; i < 30; ++i) {
            bool hasExtrema = findExtrema_(minVm, maxVm, minP, maxP, a, b, Tmin);
            if (hasExtrema)
                break;
            Tmin /= 2;
        };

        Evaluation T = Tmin;

        // Newton's method: Start at minimum temperature and minimize
        // the "gap" between the extrema of the EOS
        unsigned iMax = 100;
        for (unsigned i = 0; i < iMax; ++i) {
            // calculate function we would like to minimize
            Evaluation f = maxVm - minVm;

            // check if we're converged
            if (f < 1e-10 || (i == iMax - 1 && f < 1e-8)) {
                Tcrit = T;
                pcrit = (minP + maxP)/2;
                Vcrit = (maxVm + minVm)/2;
                return;
            }

            // backward differences. Forward differences are not
            // robust, since the isotherm could be critical if an
            // epsilon was added to the temperature. (this is case
            // rarely happens, though)
            const Scalar eps = - 1e-11;
            bool hasExtrema OPM_OPTIM_UNUSED = findExtrema_(minVm, maxVm, minP, maxP, a, b, T + eps);
            assert(hasExtrema);
            assert(std::isfinite(Opm::scalarValue(maxVm)));
            Evaluation fStar = maxVm - minVm;

            // derivative of the difference between the maximum's
            // molar volume and the minimum's molar volume regarding
            // temperature
            Evaluation fPrime = (fStar - f)/eps;
            if (std::abs(Opm::scalarValue(fPrime)) < 1e-40) {
                Tcrit = T;
                pcrit = (minP + maxP)/2;
                Vcrit = (maxVm + minVm)/2;
                return;
            }

            // update value for the current iteration
            Evaluation delta = f/fPrime;
            assert(std::isfinite(Opm::scalarValue(delta)));
            if (delta > 0)
                delta = -10;

            // line search (we have to make sure that both extrema
            // still exist after the update)
            for (unsigned j = 0; ; ++j) {
                if (j >= 20) {
                    if (f < 1e-8) {
                        Tcrit = T;
                        pcrit = (minP + maxP)/2;
                        Vcrit = (maxVm + minVm)/2;
                        return;
                    }

                    std::ostringstream oss;
                    oss << "Could not determine the critical point for a=" << a << ", b=" << b;
                    throw NumericalIssue(oss.str());
                }

                if (findExtrema_(minVm, maxVm, minP, maxP, a, b, T - delta)) {
                    // if the isotherm for T - delta exhibits two
                    // extrema the update is finished
                    T -= delta;
                    break;
                }
                else
                    delta /= 2;
            };
        }

        assert(false);
    }

    // find the two molar volumes where the EOS exhibits extrema and
    // which are larger than the covolume of the phase
    template <class Evaluation>
    static bool findExtrema_(Evaluation& Vmin,
                             Evaluation& Vmax,
                             Evaluation& /*pMin*/,
                             Evaluation& /*pMax*/,
                             const Evaluation& a,
                             const Evaluation& b,
                             const Evaluation& T)
    {
        Scalar u = 2;
        Scalar w = -1;

        const Evaluation& RT = R*T;

        // calculate coefficients of the 4th order polynominal in
        // monomial basis
        const Evaluation& a1 = RT;
        const Evaluation& a2 = 2*RT*u*b - 2*a;
        const Evaluation& a3 = 2*RT*w*b*b + RT*u*u*b*b  + 4*a*b - u*a*b;
        const Evaluation& a4 = 2*RT*u*w*b*b*b + 2*u*a*b*b - 2*a*b*b;
        const Evaluation& a5 = RT*w*w*b*b*b*b - u*a*b*b*b;

        assert(std::isfinite(Opm::scalarValue(a1)));
        assert(std::isfinite(Opm::scalarValue(a2)));
        assert(std::isfinite(Opm::scalarValue(a3)));
        assert(std::isfinite(Opm::scalarValue(a4)));
        assert(std::isfinite(Opm::scalarValue(a5)));

        // Newton method to find first root

        // if the values which we got on Vmin and Vmax are usefull, we
        // will reuse them as initial value, else we will start 10%
        // above the covolume
        Evaluation V = b*1.1;
        Evaluation delta = 1.0;
        for (unsigned i = 0; std::abs(Opm::scalarValue(delta)) > 1e-12; ++i) {
            const Evaluation& f = a5 + V*(a4 + V*(a3 + V*(a2 + V*a1)));
            const Evaluation& fPrime = a4 + V*(2*a3 + V*(3*a2 + V*4*a1));

            if (std::abs(Opm::scalarValue(fPrime)) < 1e-20) {
                // give up if the derivative is zero
                return false;
            }


            delta = f/fPrime;
            V -= delta;

            if (i > 200) {
                // give up after 200 iterations...
                return false;
            }
        }
        assert(std::isfinite(Opm::scalarValue(V)));

        // polynomial division
        Evaluation b1 = a1;
        Evaluation b2 = a2 + V*b1;
        Evaluation b3 = a3 + V*b2;
        Evaluation b4 = a4 + V*b3;

        // invert resulting cubic polynomial analytically
        Evaluation allV[4];
        allV[0] = V;
        int numSol = 1 + Opm::invertCubicPolynomial<Evaluation>(allV + 1, b1, b2, b3, b4);

        // sort all roots of the derivative
        std::sort(allV + 0, allV + numSol);

        // check whether the result is physical
        if (numSol != 4 || allV[numSol - 2] < b) {
            // the second largest extremum is smaller than the phase's
            // covolume which is physically impossible
            return false;
        }


        // it seems that everything is okay...
        Vmin = allV[numSol - 2];
        Vmax = allV[numSol - 1];
        return true;
    }

    /*!
     * \brief The Ambrose-Walton method to estimate the vapor
     *        pressure.
     *
     * \return Vapor pressure estimate in bar
     *
     * See:
     *
     * D. Ambrose, J. Walton: "Vapor Pressures up to Their Critical
     * Temperatures of Normal Alkanes and 1-Alkanols", Pure
     * Appl. Chem., 61, 1395-1403, 1989
     */
    template <class Evaluation, class Params>
    static Evaluation ambroseWalton_(const Params& /*params*/, const Evaluation& T)
    {
        typedef typename Params::Component Component;

        const Evaluation& Tr = T / Component::criticalTemperature();
        const Evaluation& tau = 1 - Tr;
        const Evaluation& omega = Component::acentricFactor();

        const Evaluation& f0 = (tau*(-5.97616 + Opm::sqrt(tau)*(1.29874 - tau*0.60394)) - 1.06841*Opm::pow(tau, 5))/Tr;
        const Evaluation& f1 = (tau*(-5.03365 + Opm::sqrt(tau)*(1.11505 - tau*5.41217)) - 7.46628*Opm::pow(tau, 5))/Tr;
        const Evaluation& f2 = (tau*(-0.64771 + Opm::sqrt(tau)*(2.41539 - tau*4.26979)) + 3.25259*Opm::pow(tau, 5))/Tr;

        return Component::criticalPressure()*std::exp(f0 + omega * (f1 + omega*f2));
    }

    /*!
     * \brief Returns the difference between the liquid and the gas phase
     *        fugacities in [bar]
     *
     * \param params Parameters
     * \param T Temperature [K]
     * \param p Pressure [bar]
     * \param VmLiquid Molar volume of the liquid phase [cm^3/mol]
     * \param VmGas Molar volume of the gas phase [cm^3/mol]
     */
    template <class Evaluation, class Params>
    static Evaluation fugacityDifference_(const Params& params,
                                          const Evaluation& T,
                                          const Evaluation& p,
                                          const Evaluation& VmLiquid,
                                          const Evaluation& VmGas)
    { return fugacity(params, T, p, VmLiquid) - fugacity(params, T, p, VmGas); }

/*
    static UniformTabulated2DFunction<Scalar> criticalTemperature_;
    static UniformTabulated2DFunction<Scalar> criticalPressure_;
    static UniformTabulated2DFunction<Scalar> criticalMolarVolume_;
*/
};

template <class Scalar>
const Scalar PengRobinson<Scalar>::R = Opm::Constants<Scalar>::R;

/*
template <class Scalar>
UniformTabulated2DFunction<Scalar> PengRobinson<Scalar>::criticalTemperature_;

template <class Scalar>
UniformTabulated2DFunction<Scalar> PengRobinson<Scalar>::criticalPressure_;

template <class Scalar>
UniformTabulated2DFunction<Scalar> PengRobinson<Scalar>::criticalMolarVolume_;
*/

} // namespace Opm

#endif
