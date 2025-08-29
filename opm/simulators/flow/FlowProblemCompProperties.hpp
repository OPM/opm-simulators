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
 * \copydoc Opm::FlowBaseProblemComp
 */
#ifndef OPM_FLOW_PROBLEM_COMP_PROPERTIES_HPP
#define OPM_FLOW_PROBLEM_COMP_PROPERTIES_HPP


#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/FlowBaseProblemProperties.hpp>

#include <tuple>

namespace Opm {
template <class TypeTag>
class FlowProblemComp;
}

namespace Opm::Properties {

namespace TTag {

struct FlowBaseProblemComp {
    using InheritsFrom = std::tuple<FlowBaseProblem>;
};

}
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FlowBaseProblemComp>
{ using type = FlowProblemComp<TypeTag>; };

template<class TypeTag>
struct TracerModel<TypeTag, TTag::FlowBaseProblemComp> {
    using type = ::Opm::TracerModel<TypeTag>;
};

// Set the material law for fluid fluxes
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::FlowBaseProblemComp>
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

// Enable diffusion
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

} // namespace Opm::Properties
#endif // OPM_FLOW_PROBLEM_COMP_PROPERTIES_HPP
