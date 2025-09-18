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

#ifndef OPM_BLACKOILMODEL_PROPERTIES_HEADER_INCLUDED
#define OPM_BLACKOILMODEL_PROPERTIES_HEADER_INCLUDED

#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <tuple>

namespace Opm {
template<class TypeTag> class BlackoilAquiferModel;
template<class TypeTag> class BlackoilWellModel;
}

namespace Opm::Properties {

namespace TTag {

struct FlowIstlSolver;
struct FlowBaseProblemBlackoil;  
struct FlowProblem { using InheritsFrom = std::tuple<FlowBaseProblemBlackoil, BlackOilModel>; };

}

// default in flow is to formulate the equations in surface volumes
template<class TypeTag>
struct BlackoilConserveSurfaceVolume<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = true; };

template<class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct AquiferModel<TypeTag, TTag::FlowProblem>
{ using type = BlackoilAquiferModel<TypeTag>; };

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableTemperature<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableFoam<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableBrine<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableSaltPrecipitation<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableBioeffects<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = true; };

template<class TypeTag>
struct WellModel<TypeTag, TTag::FlowProblem>
{ using type = BlackoilWellModel<TypeTag>; };

template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::FlowProblem>
{ using type = TTag::FlowIstlSolver; };

template<class TypeTag>
struct EnableDebuggingChecks<TypeTag, TTag::FlowProblem>
{ static constexpr bool value = false; };

} // namespace Opm::Properties

#endif // OPM_BLACKOILMODEL_PROPERTIES_HEADER_INCLUDED
