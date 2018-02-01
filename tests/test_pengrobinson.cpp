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
 * \brief This is test for the SPE5 fluid system (which uses the
 *        Peng-Robinson EOS) and the NCP flash solver.
 */
#include "config.h"

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/Spe5FluidSystem.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

#include <dune/common/parallel/mpihelper.hh>

template <class FluidSystem, class FluidState>
void createSurfaceGasFluidSystem(FluidState& gasFluidState)
{
    static const int gasPhaseIdx = FluidSystem::gasPhaseIdx;

    // temperature
    gasFluidState.setTemperature(273.15 + 20);

    // gas pressure
    gasFluidState.setPressure(gasPhaseIdx, 1e5);

    // gas saturation
    gasFluidState.setSaturation(gasPhaseIdx, 1.0);

    //  gas composition: mostly methane, a bit of propane
    gasFluidState.setMoleFraction(gasPhaseIdx, FluidSystem::H2OIdx, 0.0);
    gasFluidState.setMoleFraction(gasPhaseIdx, FluidSystem::C1Idx, 0.94);
    gasFluidState.setMoleFraction(gasPhaseIdx, FluidSystem::C3Idx, 0.06);
    gasFluidState.setMoleFraction(gasPhaseIdx, FluidSystem::C6Idx, 0.00);
    gasFluidState.setMoleFraction(gasPhaseIdx, FluidSystem::C10Idx, 0.00);
    gasFluidState.setMoleFraction(gasPhaseIdx, FluidSystem::C15Idx, 0.00);
    gasFluidState.setMoleFraction(gasPhaseIdx, FluidSystem::C20Idx, 0.00);

    // gas density
    typename FluidSystem::template ParameterCache<typename FluidState::Scalar> paramCache;
    paramCache.updatePhase(gasFluidState, gasPhaseIdx);
    gasFluidState.setDensity(gasPhaseIdx,
                             FluidSystem::density(gasFluidState, paramCache, gasPhaseIdx));
}

template <class Scalar, class FluidSystem, class FluidState>
Scalar computeSumxg(FluidState& resultFluidState,
                    const FluidState& prestineFluidState,
                    const FluidState& gasFluidState,
                    Scalar additionalGas)
{
    static const int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static const int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static const int numComponents = FluidSystem::numComponents;

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Opm::NcpFlash<Scalar, FluidSystem> Flash;

    resultFluidState.assign(prestineFluidState);

    // add a bit of additional gas components
    ComponentVector totalMolarities;
    for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++ compIdx)
        totalMolarities =
            prestineFluidState.molarity(oilPhaseIdx, compIdx)
            + additionalGas*gasFluidState.moleFraction(gasPhaseIdx, compIdx);

    // "flash" the modified fluid state
    typename FluidSystem::ParameterCache paramCache;
    Flash::solve(resultFluidState, totalMolarities);

    Scalar sumxg = 0;
    for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx)
        sumxg += resultFluidState.moleFraction(gasPhaseIdx, compIdx);

    return sumxg;
}

template <class Scalar, class FluidSystem, class FluidState>
void makeOilSaturated(FluidState& fluidState, const FluidState& gasFluidState)
{
    static const int gasPhaseIdx = FluidSystem::gasPhaseIdx;

    FluidState prestineFluidState;
    prestineFluidState.assign(fluidState);

    Scalar sumxg = 0;
    for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx)
        sumxg += fluidState.moleFraction(gasPhaseIdx, compIdx);

    // Newton method
    Scalar tol = 1e-8;
    Scalar additionalGas = 0; // [mol]
    for (int i = 0; std::abs(sumxg - 1) > tol; ++i) {
        if (i > 50)
            throw std::runtime_error("Newton method did not converge after 50 iterations");

        Scalar eps = std::max(1e-8, additionalGas*1e-8);

        Scalar f = 1 - computeSumxg<Scalar, FluidSystem>(prestineFluidState,
                                                         fluidState,
                                                         gasFluidState,
                                                         additionalGas);
        Scalar fStar = 1 - computeSumxg<Scalar, FluidSystem>(prestineFluidState,
                                                             fluidState,
                                                             gasFluidState,
                                                             additionalGas + eps);
        Scalar fPrime = (fStar - f)/eps;

        additionalGas -= f/fPrime;
    };
}

template <class FluidSystem, class FluidState>
void guessInitial(FluidState& fluidState, unsigned phaseIdx)
{
    if (phaseIdx == FluidSystem::gasPhaseIdx) {
        fluidState.setMoleFraction(phaseIdx, FluidSystem::H2OIdx, 0.0);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C1Idx, 0.74785);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C3Idx, 0.0121364);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C6Idx, 0.00606028);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C10Idx, 0.00268136);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C15Idx, 0.000204256);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C20Idx, 8.78291e-06);
    }
    else if (phaseIdx == FluidSystem::oilPhaseIdx) {
        fluidState.setMoleFraction(phaseIdx, FluidSystem::H2OIdx, 0.0);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C1Idx, 0.50);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C3Idx, 0.03);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C6Idx, 0.07);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C10Idx, 0.20);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C15Idx, 0.15);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::C20Idx, 0.05);
    }
    else {
        assert(phaseIdx == FluidSystem::waterPhaseIdx);
    }
}

template <class Scalar, class FluidSystem, class FluidState>
Scalar bringOilToSurface(FluidState& surfaceFluidState, Scalar alpha, const FluidState& reservoirFluidState, bool guessInitial)
{
    enum {
        numPhases = FluidSystem::numPhases,
        waterPhaseIdx = FluidSystem::waterPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
        oilPhaseIdx = FluidSystem::oilPhaseIdx,

        numComponents = FluidSystem::numComponents
    };

    typedef Opm::NcpFlash<Scalar, FluidSystem> Flash;
    typedef Opm::ThreePhaseMaterialTraits<Scalar, waterPhaseIdx, oilPhaseIdx, gasPhaseIdx> MaterialTraits;
    typedef Opm::LinearMaterial<MaterialTraits> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    const Scalar refPressure = 1.0135e5; // [Pa]

    // set the parameters for the capillary pressure law
    MaterialLawParams matParams;
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        matParams.setPcMinSat(phaseIdx, 0.0);
        matParams.setPcMaxSat(phaseIdx, 0.0);
    }
    matParams.finalize();

    // retieve the global volumetric component molarities
    surfaceFluidState.setTemperature(273.15 + 20);

    ComponentVector molarities;
    for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
        molarities[compIdx] = reservoirFluidState.molarity(oilPhaseIdx, compIdx);

    if (guessInitial) {
        // we start at a fluid state with reservoir oil.
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
                surfaceFluidState.setMoleFraction(phaseIdx,
                                                  compIdx,
                                                  reservoirFluidState.moleFraction(phaseIdx, compIdx));
            }
            surfaceFluidState.setDensity(phaseIdx, reservoirFluidState.density(phaseIdx));
            surfaceFluidState.setPressure(phaseIdx, reservoirFluidState.pressure(phaseIdx));
            surfaceFluidState.setSaturation(phaseIdx, 0.0);
        }
        surfaceFluidState.setSaturation(oilPhaseIdx, 1.0);
        surfaceFluidState.setSaturation(gasPhaseIdx, 1.0 - surfaceFluidState.saturation(oilPhaseIdx));
    }

    typename FluidSystem::template ParameterCache<Scalar> paramCache;
    paramCache.updateAll(surfaceFluidState);

    // increase volume until we are at surface pressure. use the
    // newton method for this
    ComponentVector tmpMolarities;
    for (int i = 0;; ++i) {
        if (i >= 20)
            throw Opm::NumericalIssue("Newton method did not converge after 20 iterations");

        // calculate the deviation from the standard pressure
        tmpMolarities = molarities;
        tmpMolarities /= alpha;
        Flash::template solve<MaterialLaw>(surfaceFluidState, matParams, paramCache, tmpMolarities);
        Scalar f = surfaceFluidState.pressure(gasPhaseIdx) - refPressure;

        // calculate the derivative of the deviation from the standard
        // pressure
        Scalar eps = alpha*1e-10;
        tmpMolarities = molarities;
        tmpMolarities /= alpha + eps;
        Flash::template solve<MaterialLaw>(surfaceFluidState, matParams, paramCache, tmpMolarities);
        Scalar fStar = surfaceFluidState.pressure(gasPhaseIdx) - refPressure;
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
    Flash::template solve<MaterialLaw>(surfaceFluidState, matParams, paramCache, tmpMolarities);
    return alpha;
}

template <class RawTable>
void printResult(const RawTable& rawTable,
                 const std::string& fieldName,
                 size_t firstIdx,
                 size_t secondIdx,
                 double hiresThres)
{
    std::cout << "std::vector<std::pair<Scalar, Scalar> > "<<fieldName<<" = {\n";

    size_t sampleIdx = 0;
    size_t numSamples = 20;
    size_t numRawHires = 0;
    for (; rawTable[numRawHires][firstIdx] > hiresThres; ++numRawHires)
    {}

    for (; sampleIdx < numSamples; ++sampleIdx) {
        size_t rawIdx = sampleIdx*numRawHires/numSamples;
        std::cout << "{ " << rawTable[rawIdx][firstIdx] << ", "
                  << rawTable[rawIdx][secondIdx] << " }"
                  << ",\n";
    }

    numSamples = 15;
    for (sampleIdx = 0; sampleIdx < numSamples; ++sampleIdx) {
        size_t rawIdx = sampleIdx*(rawTable.size() - numRawHires)/numSamples + numRawHires;
        std::cout << "{ " << rawTable[rawIdx][firstIdx] << ", "
                  << rawTable[rawIdx][secondIdx] << " }";
        if (sampleIdx < numSamples - 1)
            std::cout << ",\n";
        else
            std::cout << "\n";
    }

    std::cout << "};\n";
}

template <class Scalar>
inline void testAll()
{
    typedef Opm::FluidSystems::Spe5<Scalar> FluidSystem;

    enum {
        numPhases = FluidSystem::numPhases,
        waterPhaseIdx = FluidSystem::waterPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
        oilPhaseIdx = FluidSystem::oilPhaseIdx,

        numComponents = FluidSystem::numComponents,
        H2OIdx = FluidSystem::H2OIdx,
        C1Idx = FluidSystem::C1Idx,
        C3Idx = FluidSystem::C3Idx,
        C6Idx = FluidSystem::C6Idx,
        C10Idx = FluidSystem::C10Idx,
        C15Idx = FluidSystem::C15Idx,
        C20Idx = FluidSystem::C20Idx
    };

    typedef Opm::NcpFlash<Scalar, FluidSystem> Flash;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Opm::CompositionalFluidState<Scalar, FluidSystem> FluidState;

    typedef Opm::ThreePhaseMaterialTraits<Scalar, waterPhaseIdx, oilPhaseIdx, gasPhaseIdx> MaterialTraits;
    typedef Opm::LinearMaterial<MaterialTraits> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    typedef typename FluidSystem::template ParameterCache<Scalar> ParameterCache;

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
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        matParams.setPcMinSat(phaseIdx, 0.0);
        matParams.setPcMaxSat(phaseIdx, 0.0);
    }
    matParams.finalize();

    ////////////
    // Create a fluid state
    ////////////
    FluidState gasFluidState;
    createSurfaceGasFluidSystem<FluidSystem>(gasFluidState);

    FluidState fluidState;
    ParameterCache paramCache;

    // temperature
    fluidState.setTemperature(T);

    // oil pressure
    fluidState.setPressure(oilPhaseIdx, 4000 * 6894.7573); // 4000 PSI

    // oil saturation
    fluidState.setSaturation(oilPhaseIdx, 1.0);
    fluidState.setSaturation(gasPhaseIdx, 1.0 - fluidState.saturation(oilPhaseIdx));

    // oil composition: SPE-5 reservoir oil
    fluidState.setMoleFraction(oilPhaseIdx, H2OIdx, 0.0);
    fluidState.setMoleFraction(oilPhaseIdx, C1Idx, 0.50);
    fluidState.setMoleFraction(oilPhaseIdx, C3Idx, 0.03);
    fluidState.setMoleFraction(oilPhaseIdx, C6Idx, 0.07);
    fluidState.setMoleFraction(oilPhaseIdx, C10Idx, 0.20);
    fluidState.setMoleFraction(oilPhaseIdx, C15Idx, 0.15);
    fluidState.setMoleFraction(oilPhaseIdx, C20Idx, 0.05);

    //makeOilSaturated<Scalar, FluidSystem>(fluidState, gasFluidState);

    // set the saturations and pressures of the other phases
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        if (phaseIdx != oilPhaseIdx) {
            fluidState.setSaturation(phaseIdx, 0.0);
            fluidState.setPressure(phaseIdx, fluidState.pressure(oilPhaseIdx));
        }

        // initial guess for the composition (needed by the ComputeFromReferencePhase
        // constraint solver. TODO: bug in ComputeFromReferencePhase?)
        guessInitial<FluidSystem>(fluidState, phaseIdx);
    }

    typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
    CFRP::solve(fluidState,
                paramCache,
                /*refPhaseIdx=*/oilPhaseIdx,
                /*setViscosity=*/false,
                /*setEnthalpy=*/false);

    ////////////
    // Calculate the total molarities of the components
    ////////////
    ComponentVector totalMolarities;
    for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
        totalMolarities[compIdx] = fluidState.saturation(oilPhaseIdx)*fluidState.molarity(oilPhaseIdx, compIdx);

    ////////////
    // Gradually increase the volume for and calculate the gas
    // formation factor, oil formation volume factor and gas formation
    // volume factor.
    ////////////

    FluidState flashFluidState, surfaceFluidState;
    flashFluidState.assign(fluidState);
    //Flash::guessInitial(flashFluidState, totalMolarities);
    Flash::template solve<MaterialLaw>(flashFluidState, matParams, paramCache, totalMolarities);

    Scalar surfaceAlpha = 1;
    surfaceAlpha = bringOilToSurface<Scalar, FluidSystem>(surfaceFluidState, surfaceAlpha, flashFluidState, /*guessInitial=*/true);
    Scalar rho_gRef = surfaceFluidState.density(gasPhaseIdx);
    Scalar rho_oRef = surfaceFluidState.density(oilPhaseIdx);

    std::vector<std::array<Scalar, 10> > resultTable;

    Scalar minAlpha = 0.98;
    Scalar maxAlpha = surfaceAlpha;

    std::cout << "alpha[-] p[Pa] S_g[-] rho_o[kg/m^3] rho_g[kg/m^3] <M_o>[kg/mol] <M_g>[kg/mol] R_s[m^3/m^3] B_g[-] B_o[-]\n";
    int n = 300;
    for (int i = 0; i < n; ++i) {
        // ratio between the original and the current volume
        Scalar alpha = minAlpha + (maxAlpha - minAlpha)*i/(n - 1);

        // increasing the volume means decreasing the molartity
        ComponentVector curTotalMolarities = totalMolarities;
        curTotalMolarities /= alpha;

        // "flash" the modified reservoir oil
        Flash::template solve<MaterialLaw>(flashFluidState, matParams, paramCache, curTotalMolarities);

        surfaceAlpha = bringOilToSurface<Scalar, FluidSystem>(surfaceFluidState,
                                                              surfaceAlpha,
                                                              flashFluidState,
                                                              /*guessInitial=*/false);
        Scalar Rs =
            surfaceFluidState.saturation(gasPhaseIdx)
            / surfaceFluidState.saturation(oilPhaseIdx);
        std::cout << alpha << " "
                  << flashFluidState.pressure(oilPhaseIdx) << " "
                  << flashFluidState.saturation(gasPhaseIdx) << " "
                  << flashFluidState.density(oilPhaseIdx) << " "
                  << flashFluidState.density(gasPhaseIdx) << " "
                  << flashFluidState.averageMolarMass(oilPhaseIdx) << " "
                  << flashFluidState.averageMolarMass(gasPhaseIdx) << " "
                  << Rs << " "
                  << rho_gRef/flashFluidState.density(gasPhaseIdx) << " "
                  << rho_oRef/flashFluidState.density(oilPhaseIdx) << " "
                  << "\n";

        std::array<Scalar, 10> tmp;
        tmp[0] = alpha;
        tmp[1] = flashFluidState.pressure(oilPhaseIdx);
        tmp[2] = flashFluidState.saturation(gasPhaseIdx);
        tmp[3] = flashFluidState.density(oilPhaseIdx);
        tmp[4] = flashFluidState.density(gasPhaseIdx);
        tmp[5] = flashFluidState.averageMolarMass(oilPhaseIdx);
        tmp[6] = flashFluidState.averageMolarMass(gasPhaseIdx);
        tmp[7] = Rs;
        tmp[8] = rho_gRef/flashFluidState.density(gasPhaseIdx);
        tmp[9] = rho_oRef/flashFluidState.density(oilPhaseIdx);

        resultTable.push_back(tmp);
    }

    std::cout << "reference density oil [kg/m^3]: " << rho_oRef << "\n";
    std::cout << "reference density gas [kg/m^3]: " << rho_gRef << "\n";

    Scalar hiresThresholdPressure = resultTable[20][1];
    printResult(resultTable,
                "Bg", /*firstIdx=*/1, /*secondIdx=*/8,
                /*hiresThreshold=*/hiresThresholdPressure);
    printResult(resultTable,
                "Bo", /*firstIdx=*/1, /*secondIdx=*/9,
                /*hiresThreshold=*/hiresThresholdPressure);
    printResult(resultTable,
                "Rs", /*firstIdx=*/1, /*secondIdx=*/7,
                /*hiresThreshold=*/hiresThresholdPressure);
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();

    // the Peng-Robinson test currently does not work with single-precision floating
    // point scalars because of precision issues. (these are caused by the fact that the
    // test uses finite differences to calculate derivatives.)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunreachable-code"
    while (0) testAll<float>();
#pragma GCC diagnostic pop

    return 0;
}
