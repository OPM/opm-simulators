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
 * \brief Contains the classes required to extend the black-oil model by energy.
 */
#ifndef OPM_BLACK_OIL_ENERGY_MODULE_GLOBAL_INDEX_HH
#define OPM_BLACK_OIL_ENERGY_MODULE_GLOBAL_INDEX_HH

#include <opm/models/blackoil/blackoilenergymodules.hh>

namespace Opm {

template <class TypeTag, bool enableEnergyV>
class BlackOilEnergyIntensiveQuantitiesGlobalIndex;

/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by energy using global indices.
 */
template <class TypeTag>
class BlackOilEnergyIntensiveQuantitiesGlobalIndex<TypeTag, true>
    : public BlackOilEnergyIntensiveQuantities<TypeTag,true>
{
    using Parent =  BlackOilEnergyIntensiveQuantities<TypeTag, true>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidEnergyLaw = GetPropType<TypeTag, Properties::SolidEnergyLaw>;
    using ThermalConductionLaw  = GetPropType<TypeTag, Properties::ThermalConductionLaw>;
    using ParamCache = typename FluidSystem::template ParameterCache<Evaluation>;
    static constexpr bool enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>();

    using Indices = GetPropType<TypeTag, Properties::Indices>;
    static constexpr unsigned temperatureIdx = Indices::temperatureIdx;
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    void updateTemperature_([[maybe_unused]] const Problem& problem,
                            const PrimaryVariables& priVars,
                            [[maybe_unused]] unsigned globalSpaceIndex,
                            unsigned timeIdx)
    {
        auto& fs = Parent::asImp_().fluidState_;
        // set temperature
        fs.setTemperature(priVars.makeEvaluation(temperatureIdx, timeIdx));
    }
    void updateEnergyQuantities_(const Problem& problem,
                                 [[maybe_unused]] const PrimaryVariables& priVars,
                                 unsigned globalSpaceIndex,
                                 unsigned timeIdx,
                                 const ParamCache& paramCache)
    {
        auto& fs = Parent::asImp_().fluidState_;

        // compute the specific enthalpy of the fluids, the specific enthalpy of the rock
        // and the thermal conductivity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const auto& h = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
            fs.setEnthalpy(phaseIdx, h);
        }

        const auto& solidEnergyLawParams = problem().solidEnergyLawParams(globalSpaceIndex, timeIdx);
        this->rockInternalEnergy_ = SolidEnergyLaw::solidInternalEnergy(solidEnergyLawParams, fs);

        const auto& thermalConductionLawParams = problem.thermalConductionLawParams(globalSpaceIndex, timeIdx);
        this->totalThermalConductivity_ = ThermalConductionLaw::thermalConductivity(thermalConductionLawParams, fs);

        // Retrieve the rock fraction from the problem
        // Usually 1 - porosity, but if pvmult is used to modify porosity
        // we will apply the same multiplier to the rock fraction
        // i.e. pvmult*(1 - porosity) and thus interpret multpv as a volume
        // multiplier. This is to avoid negative rock volume for pvmult*porosity > 1
        this->rockFraction_ = problem.rockFraction(globalSpaceIndex, timeIdx);
    }
};

template <class TypeTag>
class BlackOilEnergyIntensiveQuantitiesGlobalIndex<TypeTag, false>
    : public BlackOilEnergyIntensiveQuantities<TypeTag, false>
{
    using Parent =  BlackOilEnergyIntensiveQuantities<TypeTag, false>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr bool enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>();

public:
    void updateTemperature_([[maybe_unused]] const Problem& problem,
                            [[maybe_unused]] const PrimaryVariables& priVars,
                            [[maybe_unused]] unsigned globalSpaceIdx,
                            [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableTemperature) {
            // even if energy is conserved, the temperature can vary over the spatial
            // domain if the EnableTemperature property is set to true
            auto& fs = this->asImp_().fluidState_;
            Scalar T = problem.temperature(globalSpaceIdx, timeIdx);
            fs.setTemperature(T);
        }
    }

    void updateEnergyQuantities_([[maybe_unused]] const Problem& problem,
                                 [[maybe_unused]] const PrimaryVariables& priVars,
                                 [[maybe_unused]] unsigned globalSpaceIdx,
                                 [[maybe_unused]] unsigned timeIdx,
                                 const typename FluidSystem::template ParameterCache<Evaluation>&)
    { }
};

} // namespace Opm

#endif
