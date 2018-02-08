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
 * \brief This is a program to test the flash calculation which uses
 *        non-linear complementarity problems (NCP)
 *
 * A flash calculation determines the pressures, saturations and
 * composition of all phases given the total mass (or, as in this case
 * the total number of moles) in a given amount of pore space.
 */
#include "config.h"

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/constraintsolvers/MiscibleMultiPhaseComposition.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/constraintsolvers/ImmiscibleFlash.hpp>

#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>

#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <sstream>

template <class Scalar, class FluidState>
void checkSame(const FluidState& fsRef, const FluidState& fsFlash)
{
    enum { numPhases = FluidState::numPhases };
    enum { numComponents = FluidState::numComponents };

    Scalar tol = std::max(std::numeric_limits<Scalar>::epsilon()*1e4, 1e-6);

    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        Scalar error;

        // check the pressures
        error = 1 - fsRef.pressure(phaseIdx)/fsFlash.pressure(phaseIdx);
        if (std::abs(error) > tol) {
                std::ostringstream oss;
                oss << "pressure error phase " << phaseIdx << " is incorrect: "
                    << fsFlash.pressure(phaseIdx)  << " flash vs "
                    << fsRef.pressure(phaseIdx) << " reference"
                    << " error=" << error;
                throw std::runtime_error(oss.str());
        }

        // check the saturations
        error = fsRef.saturation(phaseIdx) - fsFlash.saturation(phaseIdx);
        if (std::abs(error) > tol) {
            std::ostringstream oss;
            oss << "saturation error phase " << phaseIdx << " is incorrect: "
                << fsFlash.saturation(phaseIdx) << " flash vs "
                << fsRef.saturation(phaseIdx) << " reference"
                << " error=" << error;
            throw std::runtime_error(oss.str());
        }

        // check the compositions
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
            error = fsRef.moleFraction(phaseIdx, compIdx) - fsFlash.moleFraction(phaseIdx, compIdx);
            if (std::abs(error) > tol) {
                std::ostringstream oss;
                oss << "composition error phase " << phaseIdx << ", component " << compIdx
                    << " is incorrect: "
                    << fsFlash.moleFraction(phaseIdx, compIdx) << " flash vs "
                    << fsRef.moleFraction(phaseIdx, compIdx) << " reference"
                    << " error=" << error;
                throw std::runtime_error(oss.str());
            }
        }
    }
}

template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void checkImmiscibleFlash(const FluidState& fsRef,
                          typename MaterialLaw::Params& matParams)
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    // calculate the total amount of stuff in the reference fluid
    // phase
    ComponentVector globalMolarities(0.0);
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            globalMolarities[compIdx] +=
                fsRef.saturation(phaseIdx)*fsRef.molarity(phaseIdx, compIdx);
        }
    }

    // initialize the fluid state for the flash calculation
    typedef Opm::ImmiscibleFlash<Scalar, FluidSystem> ImmiscibleFlash;
    FluidState fsFlash;

    fsFlash.setTemperature(fsRef.temperature(/*phaseIdx=*/0));

    // run the flash calculation
    ImmiscibleFlash::guessInitial(fsFlash, globalMolarities);
    typename FluidSystem::template ParameterCache<typename FluidState::Scalar> paramCache;
    ImmiscibleFlash::template solve<MaterialLaw>(fsFlash, matParams, paramCache, globalMolarities);

    // compare the "flashed" fluid state with the reference one
    checkSame<Scalar>(fsRef, fsFlash);
}


template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void completeReferenceFluidState(FluidState& fs,
                                 typename MaterialLaw::Params& matParams,
                                 unsigned refPhaseIdx)
{
    enum { numPhases = FluidSystem::numPhases };
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

    unsigned otherPhaseIdx = 1 - refPhaseIdx;

    // calculate the other saturation
    fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

    // calulate the capillary pressure
    PhaseVector pC;
    MaterialLaw::capillaryPressures(pC, matParams, fs);
    fs.setPressure(otherPhaseIdx,
                   fs.pressure(refPhaseIdx)
                   + (pC[otherPhaseIdx] - pC[refPhaseIdx]));

    // set all phase densities
    typename FluidSystem::template ParameterCache<typename FluidState::Scalar> paramCache;
    paramCache.updateAll(fs);
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
        fs.setDensity(phaseIdx, rho);
    }
}

template <class Scalar>
inline void testAll()
{
    typedef Opm::FluidSystems::H2ON2<Scalar> FluidSystem;
    typedef Opm::ImmiscibleFluidState<Scalar, FluidSystem> ImmiscibleFluidState;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    enum { H2OIdx = FluidSystem::H2OIdx };
    enum { N2Idx = FluidSystem::N2Idx };

    typedef Opm::TwoPhaseMaterialTraits<Scalar, liquidPhaseIdx, gasPhaseIdx> MaterialLawTraits;
    typedef Opm::RegularizedBrooksCorey<MaterialLawTraits> EffMaterialLaw;
    typedef Opm::EffToAbsLaw<EffMaterialLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    std::cout << "---- using " << Dune::className<Scalar>() << " as scalar ----\n";
    Scalar T = 273.15 + 25;

    // initialize the tables of the fluid system
    Scalar Tmin = T - 1.0;
    Scalar Tmax = T + 1.0;
    unsigned nT = 3;

    Scalar pmin = 0.0;
    Scalar pmax = 1.25 * 2e6;
    unsigned np = 100;

    FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);

    // set the parameters for the capillary pressure law
    MaterialLawParams matParams;
    matParams.setResidualSaturation(MaterialLaw::wettingPhaseIdx, 0.0);
    matParams.setResidualSaturation(MaterialLaw::nonWettingPhaseIdx, 0.0);
    matParams.setEntryPressure(0);
    matParams.setLambda(2.0);
    matParams.finalize();

    ImmiscibleFluidState fsRef;

    // create an fluid state which is consistent

    // set the fluid temperatures
    fsRef.setTemperature(T);

    ////////////////
    // only liquid
    ////////////////
    std::cout << "testing single-phase liquid\n";

    // set liquid saturation and pressure
    fsRef.setSaturation(liquidPhaseIdx, 1.0);
    fsRef.setPressure(liquidPhaseIdx, 1e6);

    // set the remaining parameters of the reference fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams, liquidPhaseIdx);

    // check the flash calculation
    checkImmiscibleFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // only gas
    ////////////////
    std::cout << "testing single-phase gas\n";

    // set gas saturation and pressure
    fsRef.setSaturation(gasPhaseIdx, 1.0);
    fsRef.setPressure(gasPhaseIdx, 1e6);

    // set the remaining parameters of the reference fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams, gasPhaseIdx);

    // check the flash calculation
    checkImmiscibleFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // both phases
    ////////////////
    std::cout << "testing two-phase\n";

    // set liquid saturation and pressure
    fsRef.setSaturation(liquidPhaseIdx, 0.5);
    fsRef.setPressure(liquidPhaseIdx, 1e6);

    // set the remaining parameters of the reference fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams, liquidPhaseIdx);

    // check the flash calculation
    checkImmiscibleFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // with capillary pressure
    ////////////////
    std::cout << "testing two-phase with capillary pressure\n";

    MaterialLawParams matParams2;
    matParams2.setResidualSaturation(MaterialLaw::wettingPhaseIdx, 0.0);
    matParams2.setResidualSaturation(MaterialLaw::nonWettingPhaseIdx, 0.0);
    matParams2.setEntryPressure(1e3);
    matParams2.setLambda(2.0);
    matParams2.finalize();

    // set liquid saturation
    fsRef.setSaturation(liquidPhaseIdx, 0.5);

    // set pressure of the liquid phase
    fsRef.setPressure(liquidPhaseIdx, 1e6);

    // set the remaining parameters of the reference fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams2, liquidPhaseIdx);

    // check the flash calculation
    checkImmiscibleFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams2);
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    testAll<float>();

    return 0;
}
