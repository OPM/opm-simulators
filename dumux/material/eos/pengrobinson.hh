// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
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
#ifndef DUMUX_PENG_ROBINSON_HH
#define DUMUX_PENG_ROBINSON_HH

#include <dumux/material/idealgas.hh>
#include <dumux/common/math.hh>

#include <iostream>

namespace Dumux
{

/*!
 *
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
    static constexpr Scalar R = Dumux::Constants<Scalar>::R;

    PengRobinson()
    { }

public:
    /*!
     * \brief Predicts the vapor pressure for the temperature given in
     *        setTP().
     *
     * Initially, the vapor pressure is roughly estimated by using the
     * Ambrose-Walton method, then the Newton method is used to make
     * difference between the gas and liquid phase fugacity zero.
     */
    template <class Params>
    static Scalar computeVaporPressure(const Params &params, Scalar T)
    {
        typedef typename Params::Component Component;
        if (T >= Component::criticalTemperature())
            return Component::criticalPressure();

        // initial guess of the vapor pressure
        Scalar Vm[3];
        const Scalar eps = Component::criticalPressure()*1e-10;

        // use the Ambrose-Walton method to get an initial guess of
        // the vapor pressure
        Scalar pVap = ambroseWalton_(params, T);

        // Newton-Raphson method
        for (int i = 0; i < 5; ++i) {
            // calculate the molar densities
            int numSol = molarVolumes(Vm, params, T, pVap);
            assert(numSol == 3);

            Scalar f = fugacityDifference_(params, T, pVap, Vm[0], Vm[2]);
            Scalar df_dp =
                fugacityDifference_(params, T, pVap  + eps, Vm[0], Vm[2])
                -
                fugacityDifference_(params, T, pVap - eps, Vm[0], Vm[2]);
            df_dp /= 2*eps;

            Scalar delta = f/df_dp;
            pVap = pVap - delta;

            if (std::abs(delta/pVap) < 1e-10)
                break;
        }

        return pVap;
    }

    /*!
     * \brief Computes molar volumes where the Peng-Robinson EOS is
     *        true.
     */
    template <class Params>
    static Scalar computeMolarVolume(const Params &params,
                                     int phaseIdx,
                                     bool gasPhase)
    {
        Valgrind::CheckDefined(params.temperature(phaseIdx));
        Valgrind::CheckDefined(params.pressure(phaseIdx));

        Scalar Vm;
        Valgrind::SetUndefined(Vm);

        Scalar T = params.temperature(phaseIdx);
        Scalar p = params.pressure(phaseIdx);

        Scalar RT= R*T;
        Scalar Astar = params.a(phaseIdx)*p/(RT*RT);
        Scalar Bstar = params.b(phaseIdx)*p/RT;

        Scalar a1 = 1.0;
        Scalar a2 = - (1 - Bstar);
        Scalar a3 = Astar - Bstar*(3*Bstar + 2);
        Scalar a4 = Bstar*(- Astar + Bstar*(1 + Bstar));

        // ignore the first two results if the smallest
        // compressibility factor is <= 0.0. (this means that if we
        // would get negative molar volumes for the liquid phase, we
        // consider the liquid phase non-existant.)
        Scalar Z[3] = {0.0,0.0,0.0};
        int numSol = invertCubicPolynomial(Z, a1, a2, a3, a4);
        if (numSol == 3) {
            // the EOS has three intersections with the pressure,
            // i.e. the molar volume of gas is the largest one and the
            // molar volume of liquid is the smallest one
            if (gasPhase)
                Vm = Z[2]*RT/p;
            else
                Vm = Z[0]*RT/p;
        }
        else if (numSol == 1) {
            // the EOS only has one intersection with the pressure,
            // for the other phase, we take the extremum of the EOS
            // with the largest distance from the intersection.
            Scalar VmCubic = Z[0]*RT/p;

            // find the extrema (if they are present)
            Scalar Vmin, Vmax;
            bool hasExtrema;
            hasExtrema = findExtrema_(Vmin, Vmax, params, phaseIdx);

            if (!hasExtrema) {
                // if the EOS does not exhibit any extrema, the fluid
                // is critical...
                Vm = VmCubic;
                handleCriticalFluid_(Vm, params, phaseIdx, gasPhase);
            }
            else {
                if (gasPhase)
                    Vm = std::max(Vmax, VmCubic);
                else
                    Vm = std::min(Vmin, VmCubic);
            }
        }

        Valgrind::CheckDefined(Vm);
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
    template <class Params>
    static Scalar computeFugacityCoeff(const Params &params)
    {
        Scalar T = params.temperature();
        Scalar p = params.pressure();
        Scalar Vm = params.molarVolume();

        Scalar RT = R*T;
        Scalar Z = p*Vm/RT;
        Scalar Bstar = p*params.b() / RT;

        Scalar tmp =
            (Vm + params.b()*(1 + std::sqrt(2))) /
            (Vm + params.b()*(1 - std::sqrt(2)));
        Scalar expo = - params.a()/(RT * 2 * params.b() * std::sqrt(2));
        Scalar fugCoeff =
            std::exp(Z - 1) / (Z - Bstar) *
            std::pow(tmp, expo);

        return fugCoeff;
    };

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
    template <class Params>
    static Scalar computeFugacity(const Params &params)
    { return params.pressure()*computeFugacityCoeff(params); }

protected:
    template <class Params>
    static void handleCriticalFluid_(Scalar &Vm,
                                     const Params &params,
                                     int phaseIdx,
                                     bool gasPhase)
    {
        Scalar Vcrit;
        findCriticalMolarVolume_(Vcrit,
                                 params,
                                 phaseIdx,
                                 Vm,
                                 gasPhase);
        if (gasPhase)
            Vm = std::max(Vm, Vcrit);
        else
            Vm = std::min(Vm, Vcrit);
    }

    template <class Params>
    static void findCriticalMolarVolume_(Scalar &Vcrit,
                                         const Params &params,
                                         int phaseIdx,
                                         Scalar Vcubic,
                                         bool gasPhase)
    {
        Params tmpParams(params);

        // use an isotherm where the EOS should exhibit two maxima
        Scalar Tmin = 150; // [K]

        Scalar minVm;
        Scalar maxVm;
        Scalar T = Tmin;
        for (int i = 0; ; ++i) {
            tmpParams.setTemperature(phaseIdx, T);
            tmpParams.updateEosParams(phaseIdx);

            if (findExtrema_(minVm, maxVm, tmpParams, phaseIdx, /*hintSet=*/false))
                break;

            T = (T + params.temperature(phaseIdx)) / 2;
            if (i >= 3) {
                DUNE_THROW(NumericalProblem,
                           "TODO: better way to identify a non-critical temperature");
            };
        }

        // Newton's method: Start at minimum temperature and minimize
        // the "gap" between the extrema of the EOS
        for (int i = 0; i < 25; ++i) {
            // calculate function we would like to minimize
            Scalar f = maxVm - minVm;

            // check if we're converged
            if (f < 1e-8 ||
                (gasPhase && maxVm < Vcubic) ||
                (!gasPhase && minVm > Vcubic))
            {
                Vcrit = (maxVm + minVm)/2;
                return;
            }

            // backward differences. Forward differences are not
            // robust, since the isotherm could be critical if an
            // epsilon was added to the temperature. (this is case
            // rarely happens, though)
            const Scalar eps = - 1e-8;
            T = T + eps;
            tmpParams.updateEosParams(phaseIdx);
            tmpParams.setTemperature(phaseIdx, T);
            if (!findExtrema_(minVm, maxVm, tmpParams, phaseIdx, /*hintSet=*/true))
                DUNE_THROW(NumericalProblem,
                           "Error when calculating derivative of extrema "
                           "to calculate the critical point of phase "
                           << phaseIdx);
            T = T - eps;
            Scalar fStar = maxVm - minVm;

            // derivative of the difference between the maximum's
            // molar volume and the minimum's molar volume regarding
            // temperature
            Scalar fPrime = (fStar - f)/eps;

            // update value for the current iteration
            Scalar delta = f/fPrime;
            if (delta > 0)
                delta = -10;

            // line search (we have to make sure that both extrema
            // still exist after the update)
            Scalar Tstar = T - delta;
            for (int j = 0; ; ++j) {
                if (j >= 10) {
                    DUNE_THROW(NumericalProblem,
                               "Could not determine the critical point of phase "
                               << phaseIdx);
                }

                tmpParams.setTemperature(phaseIdx, Tstar);
                tmpParams.updateEosParams(phaseIdx);
                if (findExtrema_(minVm, maxVm, tmpParams, phaseIdx, /*hintSet=*/(j==0))) {
                    // if the isotherm at Tstar exhibits two extrema
                    // the update is finished
                    T = Tstar;
                    break;
                }
                else {
                    delta /= 2;
                    Tstar = T - delta;
                }
            };
        }
    };

    // find the two molar volumes where the EOS exhibits extrema and
    // which are larger than the covolume of the phase
    template <class Params>
    static bool findExtrema_(Scalar &Vmin,
                             Scalar &Vmax,
                             const Params &params,
                             int phaseIdx,
                             bool hintSet = false)
    {
        Scalar a = params.a(phaseIdx);
        Scalar b = params.b(phaseIdx);
        Scalar u = 2;
        Scalar w = -1;

        Scalar RT = R*params.temperature(phaseIdx);

        // calculate coefficients of the 4th order polynominal in
        // monomial basis
        Scalar a1 = RT;
        Scalar a2 = 2*RT*u*b - 2*a;
        Scalar a3 = 2*RT*w*b*b + RT*u*u*b*b  + 4*a*b - u*a*b;
        Scalar a4 = 2*RT*u*w*b*b*b + 2*u*a*b*b - 2*a*b*b;
        Scalar a5 = RT*w*w*b*b*b*b - u*a*b*b*b;

        // Newton method to find first root

        // if the values which we got on Vmin and Vmax are usefull, we
        // will reuse them as initial value, else we will start 10%
        // above the covolume
        Scalar V = hintSet?Vmin:b*1.1;
        Scalar delta = 1.0;
        for (int i = 0; std::abs(delta) > 1e-9; ++i) {
            Scalar f = a5 + V*(a4 + V*(a3 + V*(a2 + V*a1)));
            Scalar fPrime = a4 + V*(2*a3 + V*(3*a2 + V*4*a1));

            if (std::abs(fPrime) < 1e-20) {
                // give up if the derivative is zero
                Vmin = 0;
                Vmax = 0;
                return false;
            }


            delta = f/fPrime;
            V -= delta;

            if (i > 20) {
                // give up after 20 iterations...
                Vmin = 0;
                Vmax = 0;
                return false;
            }
        }

        // polynomial division
        Scalar b1 = a1;
        Scalar b2 = a2 + V*b1;
        Scalar b3 = a3 + V*b2;
        Scalar b4 = a4 + V*b3;

        // invert resulting cubic polynomial analytically
        Scalar allV[4];
        allV[0] = V;
        int numSol = 1 + Dumux::invertCubicPolynomial(&allV[1], b1, b2, b3, b4);

        // sort all roots of the derivative
        std::sort(allV + 0, allV + numSol);

        // check whether the result is physical
        if (allV[numSol - 2] < b) {
            // the second largest extremum is smaller than the phase's
            // covolume which is physically impossible
            Vmin = Vmax = 0;
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
    template <class Params>
    static Scalar ambroseWalton_(const Params &params, Scalar T)
    {
        typedef typename Params::Component Component;

        Scalar Tr = T / Component::criticalTemperature();
        Scalar tau = 1 - Tr;
        Scalar omega = Component::acentricFactor();

        Scalar f0 = (tau*(-5.97616 + std::sqrt(tau)*(1.29874 - tau*0.60394)) - 1.06841*std::pow(tau, 5))/Tr;
        Scalar f1 = (tau*(-5.03365 + std::sqrt(tau)*(1.11505 - tau*5.41217)) - 7.46628*std::pow(tau, 5))/Tr;
        Scalar f2 = (tau*(-0.64771 + std::sqrt(tau)*(2.41539 - tau*4.26979)) + 3.25259*std::pow(tau, 5))/Tr;

        return Component::criticalPressure()*std::exp(f0 + omega * (f1 + omega*f2));
    };

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
    template <class Params>
    static Scalar fugacityDifference_(const Params &params,
                                      Scalar T,
                                      Scalar p,
                                      Scalar VmLiquid,
                                      Scalar VmGas)
    { return fugacity(params, T, p, VmLiquid) - fugacity(params, T, p, VmGas); };
};

} // end namepace

#endif
