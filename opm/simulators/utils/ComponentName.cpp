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
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>
#include <opm/models/ptflash/flashindices.hh>

#include <flowexperimental/comp/flow_comp.hpp>

#include <cassert>

namespace Opm {

template<class FluidSystem, class Indices>
ComponentName<FluidSystem,Indices>::ComponentName()
    : names_(Indices::numEq)
{
    if constexpr (requires {
                      FluidSystem::solventComponentIndex(0u);
                      FluidSystem::canonicalToActiveCompIdx(0u);
                  }) {
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned canonicalCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
            names_[FluidSystem::canonicalToActiveCompIdx(canonicalCompIdx)]
                = FluidSystem::componentName(canonicalCompIdx);
        }
    }
    else {
        constexpr auto numNamedComponents =
            (Indices::numEq < FluidSystem::numComponents) ? Indices::numEq : FluidSystem::numComponents;

        for (unsigned compIdx = 0; compIdx < numNamedComponents; ++compIdx) {
            names_[compIdx] = FluidSystem::componentName(compIdx);
        }

        if constexpr (requires { Indices::waterEnabled; Indices::water0Idx; }) {
            if constexpr (Indices::waterEnabled
                          && (Indices::water0Idx >= 0)
                          && (Indices::water0Idx < Indices::numEq)) {
                names_[Indices::water0Idx] = "Water";
            }
        }
    }

    if constexpr (requires { Indices::enableSolvent; Indices::solventSaturationIdx; }) {
        if constexpr (Indices::enableSolvent) {
            names_[Indices::solventSaturationIdx] = "Solvent";
        }
    }

    if constexpr (requires { Indices::enableExtbo; Indices::zFractionIdx; }) {
        if constexpr (Indices::enableExtbo) {
            names_[Indices::zFractionIdx] = "ZFraction";
        }
    }

    if constexpr (requires { Indices::enablePolymer; Indices::polymerConcentrationIdx; }) {
        if constexpr (Indices::enablePolymer) {
            names_[Indices::polymerConcentrationIdx] = "Polymer";
        }
    }

    if constexpr (requires { Indices::polymerMoleWeightIdx; }) {
        if constexpr (Indices::polymerMoleWeightIdx >= 0) {
            names_[Indices::polymerMoleWeightIdx] = "MolecularWeightP";
        }
    }

    if constexpr (requires { Indices::enableFullyImplicitThermal; Indices::temperatureIdx; }) {
        if constexpr (Indices::enableFullyImplicitThermal) {
            names_[Indices::temperatureIdx] = "Energy";
        }
    }

    if constexpr (requires { Indices::numFoam; Indices::foamConcentrationIdx; }) {
        if constexpr (Indices::numFoam == 1) {
            names_[Indices::foamConcentrationIdx] = "Foam";
        }
    }

    if constexpr (requires { Indices::numBrine; Indices::saltConcentrationIdx; }) {
        if constexpr (Indices::numBrine == 1) {
            names_[Indices::saltConcentrationIdx] = "Brine";
        }
    }

    if constexpr (requires {
                      Indices::enableMICP;
                      Indices::microbialConcentrationIdx;
                      Indices::oxygenConcentrationIdx;
                      Indices::ureaConcentrationIdx;
                      Indices::biofilmVolumeFractionIdx;
                      Indices::calciteVolumeFractionIdx;
                  }) {
        if constexpr (Indices::enableMICP) {
            names_[Indices::microbialConcentrationIdx] = "Microbes";
            names_[Indices::oxygenConcentrationIdx] = "Oxygen";
            names_[Indices::ureaConcentrationIdx] = "Urea";
            names_[Indices::biofilmVolumeFractionIdx] = "Biofilm";
            names_[Indices::calciteVolumeFractionIdx] = "Calcite";
        }
    }

    if constexpr (requires { Indices::enableBiofilm; Indices::microbialConcentrationIdx; Indices::biofilmVolumeFractionIdx; }) {
        if constexpr (Indices::enableBiofilm) {
            names_[Indices::microbialConcentrationIdx] = "Microbes";
            names_[Indices::biofilmVolumeFractionIdx] = "Biofilm";
        }
    }
}


#include <opm/simulators/utils/InstantiationIndicesMacros.hpp>

INSTANTIATE_TYPE_INDICES(ComponentName, double)
template class ComponentName<GenericOilGasWaterFluidSystem<double, 2, false>,
                             FlashIndices<Properties::TTag::FlowCompProblem<2, false>, 0>>;
template class ComponentName<GenericOilGasWaterFluidSystem<double, 2, true>,
                             FlashIndices<Properties::TTag::FlowCompProblem<2, true>, 0>>;
template class ComponentName<GenericOilGasWaterFluidSystem<double, 3, false>,
                             FlashIndices<Properties::TTag::FlowCompProblem<3, false>, 0>>;
template class ComponentName<GenericOilGasWaterFluidSystem<double, 3, true>,
                             FlashIndices<Properties::TTag::FlowCompProblem<3, true>, 0>>;

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE_INDICES(ComponentName, float)
template class ComponentName<GenericOilGasWaterFluidSystem<float, 2, false>,
                             FlashIndices<Properties::TTag::FlowCompProblem<2, false>, 0>>;
template class ComponentName<GenericOilGasWaterFluidSystem<float, 2, true>,
                             FlashIndices<Properties::TTag::FlowCompProblem<2, true>, 0>>;
template class ComponentName<GenericOilGasWaterFluidSystem<float, 3, false>,
                             FlashIndices<Properties::TTag::FlowCompProblem<3, false>, 0>>;
template class ComponentName<GenericOilGasWaterFluidSystem<float, 3, true>,
                             FlashIndices<Properties::TTag::FlowCompProblem<3, true>, 0>>;
#endif

} // namespace Opm
