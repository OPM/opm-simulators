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
 * \copydoc Opm::FlowBaseProblemBlackoil
 */
#ifndef OPM_FLOW_BASE_PROBLEM_BLACKOIL_PROPERTIES_HPP
#define OPM_FLOW_BASE_PROBLEM_BLACKOIL_PROPERTIES_HPP


#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/FlowBaseProblemProperties.hpp>
#include <opm/simulators/flow/FIBlackoilModel.hpp>
#include <opm/simulators/flow/NewTranFluxModule.hpp>
#include <opm/simulators/flow/OutputBlackoilModule.hpp>
#include <opm/simulators/flow/BlackoilModel.hpp>


#include <tuple>

namespace Opm {
template <class TypeTag>
class FlowProblemBlackoil;
}

namespace Opm::Properties {

namespace TTag {

struct FlowBaseProblemBlackoil {
    using InheritsFrom = std::tuple<FlowBaseProblem>;
};

}

// Set the problem property
template<class TypeTag>
struct NonlinearSystem<TypeTag, TTag::FlowBaseProblemBlackoil>
{ using type = BlackoilModel<TypeTag>; };


template<class TypeTag>
struct Problem<TypeTag, TTag::FlowBaseProblemBlackoil>
{ using type = FlowProblemBlackoil<TypeTag>; };

template<class TypeTag>
struct Model<TypeTag, TTag::FlowBaseProblemBlackoil>
{ using type = FIBlackOilModel<TypeTag>; };

// Set the material law for fluid fluxes
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::FlowBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using Traits = ThreePhaseMaterialTraits<Scalar,
                                            /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                            /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                            /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx>;

public:
    using EclMaterialLawManager = ::Opm::EclMaterialLaw::Manager<Traits>;

    using type = typename EclMaterialLawManager::MaterialLaw;
};

// Use the "velocity module" which uses the Eclipse "NEWTRAN" transmissibilities
template<class TypeTag>
struct FluxModule<TypeTag, TTag::FlowBaseProblemBlackoil>
{ using type = NewTranFluxModule<TypeTag>; };

// Use the dummy gradient calculator in order not to do unnecessary work.
template<class TypeTag>
struct GradientCalculator<TypeTag, TTag::FlowBaseProblemBlackoil>
{ using type = DummyGradientCalculator<TypeTag>; };

} // namespace Opm::Properties

#endif // OPM_FLOW_BASE_PROBLEM_BLACKOIL_PROPERTIES_HPP
