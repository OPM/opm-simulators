/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/simulators/utils/ComponentName.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <cassert>

namespace Opm {

template<class FluidSystem, class Indices>
ComponentName<FluidSystem,Indices>::ComponentName()
    : names_(Indices::numEq)
{
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }

        const unsigned canonicalCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
        names_[Indices::canonicalToActiveComponentIndex(canonicalCompIdx)]
            = FluidSystem::componentName(canonicalCompIdx);
    }

    if constexpr (Indices::enableSolvent) {
        names_[Indices::solventSaturationIdx] = "Solvent";
    }

    if constexpr (Indices::enableExtbo) {
        names_[Indices::zFractionIdx] = "ZFraction";
    }

    if constexpr (Indices::enablePolymer) {
        names_[Indices::polymerConcentrationIdx] = "Polymer";
    }

    if constexpr (Indices::polymerMoleWeightIdx >= 0) {
        assert(Indices::enablePolymer);
        names_[Indices::polymerMoleWeightIdx] = "MolecularWeightP";
    }

    if constexpr (Indices::enableEnergy) {
        names_[Indices::temperatureIdx] = "Energy";
    }

    if constexpr (Indices::numFoam == 1) {
        names_[Indices::foamConcentrationIdx] = "Foam";
    }

    if constexpr (Indices::numBrine == 1) {
        names_[Indices::saltConcentrationIdx] = "Brine";
    }

    if constexpr (Indices::enableMICP) {
        names_[Indices::microbialConcentrationIdx] = "Microbes";
        names_[Indices::oxygenConcentrationIdx] = "Oxygen";
        names_[Indices::ureaConcentrationIdx] = "Urea";
        names_[Indices::biofilmConcentrationIdx] = "Biofilm";
        names_[Indices::calciteConcentrationIdx] = "Calcite";
    }
}

#define INSTANCE( ...) \
template class ComponentName<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>, \
                             __VA_ARGS__>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>)

// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

} // namespace Opm
