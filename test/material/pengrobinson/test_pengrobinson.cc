// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \brief This is test for the SPE5 fluid system (which uses the
 *        Peng-Robinson EOS) and the NCP flash solver.
 */
#include "config.h"

#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/ncpflash.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/fluidsystems/spe5fluidsystem.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mplinearmaterial.hh>

template <class FluidSystem, class FluidState>
void guessInitial(FluidState &fluidState,
                  int phaseIdx)
{
    if (phaseIdx == FluidSystem::gPhaseIdx) {
        fluidState.setMoleFraction(phaseIdx, FluidSystem::H2OIdx, 0.0);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C1Idx, 0.74785);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C3Idx, 0.0121364);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C6Idx, 0.00606028);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C10Idx, 0.00268136);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C15Idx, 0.000204256);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C20Idx, 8.78291e-06);
    }
    else if (phaseIdx == FluidSystem::oPhaseIdx) {
        fluidState.setMoleFraction(phaseIdx, FluidSystem::H2OIdx, 0.0);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C1Idx, 0.50);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C3Idx, 0.03);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C6Idx, 0.07);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C10Idx, 0.20);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C15Idx, 0.15); 
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C20Idx, 0.05); 
    }
    else {
        assert(phaseIdx == FluidSystem::wPhaseIdx);
    }
}

template <class Scalar, class FluidSystem, class FluidState>
Scalar bringOilToSurface(FluidState &surfaceFluidState, Scalar alpha, const FluidState &reservoirFluidState, bool guessInitial)
{
    enum {
        numPhases = FluidSystem::numPhases,
        wPhaseIdx = FluidSystem::wPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,
        oPhaseIdx = FluidSystem::oPhaseIdx,

        numComponents = FluidSystem::numComponents
    };

    typedef Dumux::NcpFlash<Scalar, FluidSystem> Flash;
    typedef Dumux::MpLinearMaterial<numPhases, Scalar> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    const Scalar refPressure = 1.0135e5; // [Pa]

    // set the parameters for the capillary pressure law
    MaterialLawParams matParams;
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        matParams.setPcMinSat(phaseIdx, 0.0);
        matParams.setPcMaxSat(phaseIdx, 0.0);
    }

    // retieve the global volumetric component molarities
    surfaceFluidState.setTemperature(273.15 + 20);

    ComponentVector molarities;
    for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
        molarities[compIdx] = reservoirFluidState.molarity(oPhaseIdx, compIdx);

    if (guessInitial) {
        // we start at a fluid state with reservoir oil.
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                surfaceFluidState.setMoleFraction(phaseIdx, 
                                                  compIdx,
                                                  reservoirFluidState.moleFraction(phaseIdx, compIdx));
            }
            surfaceFluidState.setDensity(phaseIdx, reservoirFluidState.density(phaseIdx));
            surfaceFluidState.setPressure(phaseIdx, reservoirFluidState.pressure(phaseIdx));
            surfaceFluidState.setSaturation(phaseIdx, 0.0);
        }
        surfaceFluidState.setSaturation(oPhaseIdx, 1.0);
    }

    typename FluidSystem::ParameterCache paramCache;
    paramCache.updateAll(surfaceFluidState);

    // increase volume until we are at surface pressure. use the
    // newton method for this
    ComponentVector tmpMolarities;
    for (int i = 0;; ++i) {
        if (i >= 20)
            DUNE_THROW(Dumux::NumericalProblem,
                       "Newton method did not converge after 20 iterations");

        // calculate the deviation from the standard pressure
        tmpMolarities = molarities;
        tmpMolarities /= alpha;
        Flash::template solve<MaterialLaw>(surfaceFluidState, paramCache, matParams, tmpMolarities);
        Scalar f = surfaceFluidState.pressure(gPhaseIdx) - refPressure;

        // calculate the derivative of the deviation from the standard
        // pressure
        Scalar eps = alpha*1e-10;
        tmpMolarities = molarities;
        tmpMolarities /= alpha + eps;
        Flash::template solve<MaterialLaw>(surfaceFluidState, paramCache, matParams, tmpMolarities);
        Scalar fStar = surfaceFluidState.pressure(gPhaseIdx) - refPressure;
        Scalar fPrime = (fStar - f)/eps;
        
        // newton update
        Scalar delta = f/fPrime;
        alpha -= delta;
        if (std::abs(delta) < std::abs(alpha)*1e-9) {
            break;
        }
    }
        
    // calculate the final result
    tmpMolarities = molarities;
    tmpMolarities /= alpha;
    Flash::template solve<MaterialLaw>(surfaceFluidState, paramCache, matParams, tmpMolarities);
    return alpha;
}

int main(int argc, char** argv)
{
    typedef double Scalar;
    typedef Dumux::FluidSystems::Spe5<Scalar> FluidSystem;

    enum {
        numPhases = FluidSystem::numPhases,
        wPhaseIdx = FluidSystem::wPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,
        oPhaseIdx = FluidSystem::oPhaseIdx,

        numComponents = FluidSystem::numComponents,
        H2OIdx = FluidSystem::H2OIdx,
        C1Idx = FluidSystem::C1Idx,
        C3Idx = FluidSystem::C3Idx,
        C6Idx = FluidSystem::C6Idx,
        C10Idx = FluidSystem::C10Idx,
        C15Idx = FluidSystem::C15Idx,
        C20Idx = FluidSystem::C20Idx
    };

    typedef Dumux::NcpFlash<Scalar, FluidSystem> Flash;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

    typedef Dumux::MpLinearMaterial<numPhases, Scalar> MaterialLaw;
    typedef MaterialLaw::Params MaterialLawParams;

    typedef FluidSystem::ParameterCache ParameterCache;

    ////////////
    // Initialize the fluid system and create the capillary pressure
    // parameters
    ////////////
    Scalar T = 273.15 + 20; // 20 deg Celsius
    FluidSystem::init(/*minTemperature=*/T - 1, 
                      /*maxTemperature=*/T + 1, 
                      /*minPressure=*/1.0e4, 
                      /*maxTemperature=*/40.0e6);

    // set the parameters for the capillary pressure law
    MaterialLawParams matParams;
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        matParams.setPcMinSat(phaseIdx, 0.0);
        matParams.setPcMaxSat(phaseIdx, 0.0);
    }

    ////////////
    // Create a fluid state
    ////////////
    FluidState fluidState;
    ParameterCache paramCache;

    // temperature
    fluidState.setTemperature(T);

    // oil pressure
    fluidState.setPressure(oPhaseIdx, 4000 * 6894.7573); // 4000 PSI

    // oil saturation
    fluidState.setSaturation(oPhaseIdx, 1.0);

    // composition: SPE-5 reservoir oil
    fluidState.setMoleFraction(oPhaseIdx, H2OIdx, 0.0);
    fluidState.setMoleFraction(oPhaseIdx, C1Idx, 0.50);
    fluidState.setMoleFraction(oPhaseIdx, C3Idx, 0.03);
    fluidState.setMoleFraction(oPhaseIdx, C6Idx, 0.07);
    fluidState.setMoleFraction(oPhaseIdx, C10Idx, 0.20);
    fluidState.setMoleFraction(oPhaseIdx, C15Idx, 0.15);
    fluidState.setMoleFraction(oPhaseIdx, C20Idx, 0.05);

    // set the fugacities in the oil phase
    paramCache.updatePhase(fluidState, oPhaseIdx);
    ComponentVector fugVec;
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        Scalar Phi =
            FluidSystem::fugacityCoefficient(fluidState, paramCache, oPhaseIdx, compIdx);
        fluidState.setFugacityCoefficient(oPhaseIdx, compIdx, Phi);
        fugVec[compIdx] = fluidState.fugacity(oPhaseIdx, compIdx);

        fluidState.setMoleFraction(gPhaseIdx, compIdx, fluidState.moleFraction(oPhaseIdx, compIdx));
    }

/*
    fluidState.setPressure(gPhaseIdx, fluidState.pressure(oPhaseIdx)); // 4000 PSI
    {
    Scalar TMin = fluidState.temperature(oPhaseIdx)/2;
    Scalar TMax = fluidState.temperature(oPhaseIdx)*2;
    int n = 1000;
    for (int i = 0; i < n; ++i) {
        Scalar T = TMin + (TMax - TMin)*i/(n - 1);

        fluidState.setTemperature(T);
        paramCache.updatePhase(fluidState, oPhaseIdx);
        paramCache.updatePhase(fluidState, gPhaseIdx);
        
        Scalar rhoO = FluidSystem::density(fluidState, paramCache, oPhaseIdx);
        Scalar rhoG = FluidSystem::density(fluidState, paramCache, gPhaseIdx);
        fluidState.setDensity(oPhaseIdx, rhoO);
        fluidState.setDensity(gPhaseIdx, rhoG);
        std::cout << T << " "
                  << fluidState.density(oPhaseIdx) << " "
                  << fluidState.density(gPhaseIdx) << " "
                  << paramCache.gasPhaseParams().a() << " "
                  << paramCache.gasPhaseParams().b() << " "
                  << "\n";
    }
    exit(0);
    }
*/

    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        if (phaseIdx != oPhaseIdx) {
            fluidState.setSaturation(phaseIdx, 0.0);
            fluidState.setPressure(phaseIdx, fluidState.pressure(oPhaseIdx));
        }

        // initial guess for the composition
        guessInitial<FluidSystem>(fluidState, phaseIdx);
    }
    
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
    CFRP::solve(fluidState,
                paramCache,
                /*refPhaseIdx=*/oPhaseIdx,
                /*setViscosity=*/false,
                /*setEnthalpy=*/false);
    
    ////////////
    // Calculate the total molarities of the components
    ////////////
    ComponentVector molarities;
    for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
        molarities[compIdx] = fluidState.saturation(oPhaseIdx)*fluidState.molarity(oPhaseIdx, compIdx);

    ////////////
    // Gradually increase the volume for and calculate the gas
    // formation factor, oil formation volume factor and gas formation
    // volume factor.
    ////////////

    FluidState flashFluidState, surfaceFluidState;
    flashFluidState.assign(fluidState);
    //Flash::guessInitial(flashFluidState, paramCache, molarities);
    Flash::solve<MaterialLaw>(flashFluidState, paramCache, matParams, molarities);

    Scalar surfaceAlpha = 1;
    surfaceAlpha = bringOilToSurface<Scalar, FluidSystem>(surfaceFluidState, surfaceAlpha, flashFluidState, /*guessInitial=*/true);
    Scalar rho_gRef = surfaceFluidState.density(gPhaseIdx);
    Scalar rho_oRef = surfaceFluidState.density(oPhaseIdx);

    std::cout << "alpha[-] p[Pa] S_g[-] rho_o[kg/m^3] rho_g[kg/m^3] <M_o>[kg/mol] <M_g>[kg/mol] R_s[m^3/m^3] B_g[-] B_o[-]\n";
    int n = 10;
    for (int i = 0; i < n; ++i) {
        Scalar minAlpha = 0.98;
        Scalar maxAlpha = 0.98 + (5 - 0.98)/100;

        // ratio between the original and the current volume
        Scalar alpha = minAlpha + (maxAlpha - minAlpha)*i/(n - 1);

        // increasing the volume means decreasing the molartity
        ComponentVector curMolarities = molarities;
        curMolarities /= alpha;

        // "flash" the modified reservoir oil
        Flash::solve<MaterialLaw>(flashFluidState, paramCache, matParams, curMolarities);

        surfaceAlpha = bringOilToSurface<Scalar, FluidSystem>(surfaceFluidState, 
                                                              surfaceAlpha, 
                                                              flashFluidState,
                                                              /*guessInitial=*/false);
        std::cout << alpha << " "
                  << flashFluidState.pressure(oPhaseIdx) << " "
                  << flashFluidState.saturation(gPhaseIdx) << " "
                  << flashFluidState.density(oPhaseIdx) << " "
                  << flashFluidState.density(gPhaseIdx) << " "
                  << flashFluidState.averageMolarMass(oPhaseIdx) << " "
                  << flashFluidState.averageMolarMass(gPhaseIdx) << " "
                  << surfaceFluidState.saturation(gPhaseIdx)*surfaceAlpha << " "
                  << flashFluidState.density(gPhaseIdx)/rho_gRef << " "
                  << flashFluidState.density(oPhaseIdx)/rho_oRef << " "
                  << "\n";
    }

    std::cout << "reference density oil [kg/m^3]: " << rho_oRef << "\n";
    std::cout << "reference density gas [kg/m^3]: " << rho_gRef << "\n";

    return 0;
}
