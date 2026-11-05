// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 Equinor ASA.

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
/*!
 * \file
 *
 * TypeTag of the black-oil TPFA flow problem without hysteresis, used by
 * \c flow for black-oil decks without HYSTER.
 *
 * Kept in a header (and not in \c flow/flow_blackoil_nohyst.cpp) so that every
 * translation unit that instantiates code with this TypeTag sees the same
 * property bindings, e.g. the GPU intensive-quantities TypeTag in
 * \c GpuFlowBlackoilNohystTypeTag.hpp.
 */
#ifndef OPM_FLOW_PROBLEM_NOHYST_TYPETAG_HPP
#define OPM_FLOW_PROBLEM_NOHYST_TYPETAG_HPP

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/BlackoilModelProperties.hpp>

#include <tuple>

namespace Opm::Properties {

namespace TTag {

struct FlowProblemNohyst
{ using InheritsFrom = std::tuple<FlowProblem>; };

}

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowProblemNohyst>
{ using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowProblemNohyst>
{ using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowProblemNohyst>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableHysteresis<TypeTag, TTag::FlowProblemNohyst>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableEndpointScaling<TypeTag, TTag::FlowProblemNohyst>
{ static constexpr bool value = true; };

template<class TypeTag>
struct AvoidElementContext<TypeTag, TTag::FlowProblemNohyst>
{ static constexpr bool value = true; };

} // namespace Opm::Properties

#endif // OPM_FLOW_PROBLEM_NOHYST_TYPETAG_HPP
