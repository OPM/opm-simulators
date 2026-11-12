/*
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
#include "config.h"

#include <flow/flow_blackoil.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicit.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>

#include <cstddef>
#include <tuple>

#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm::Properties {

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowProblemTPFA>
{ using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowProblemTPFA>
{ using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowProblemTPFA>
{ static constexpr bool value = false; };

template<class TypeTag>
struct AvoidElementContext<TypeTag, TTag::FlowProblemTPFA>
{ static constexpr bool value = true; };


namespace TTag {
// Same as FlowProblemTPFA without the convective-mixing (DRSDTCON) terms in the flux.
struct FlowProblemTPFANoConvectiveMixing
{ using InheritsFrom = std::tuple<FlowProblemTPFA>; };
}

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::FlowProblemTPFANoConvectiveMixing>
{ static constexpr bool value = false; };

} // namespace Opm::Properties

namespace {

// Convective mixing can be compiled out when no report step activates DRSDTCON.
// Decks with actions keep the general type, since an action may turn it on later.
bool deckNeedsConvectiveMixing()
{
    const auto& schedule = Opm::FlowGenericVanguard::modelParams_.eclSchedule_;
    if (!schedule) {
        return true;
    }
    for (std::size_t step = 0; step < schedule->size(); ++step) {
        const auto& state = (*schedule)[step];
        if (state.oilvap().drsdtConvective() || !state.actions().empty()) {
            return true;
        }
    }
    return false;
}

} // anonymous namespace

namespace Opm {

std::unique_ptr<FlowMain<Properties::TTag::FlowProblemTPFA>>
flowBlackoilTpfaMainInit(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    return std::make_unique<FlowMain<Properties::TTag::FlowProblemTPFA>>(
        argc, argv, outputCout, outputFiles);
}

// ----------------- Main program -----------------
int flowBlackoilTpfaMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    if (deckNeedsConvectiveMixing()) {
        FlowMain<Properties::TTag::FlowProblemTPFA>
            mainfunc {argc, argv, outputCout, outputFiles};
        return mainfunc.execute();
    }
    FlowMain<Properties::TTag::FlowProblemTPFANoConvectiveMixing>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowBlackoilTpfaMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowProblemTPFA;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStaticSelect<TypeTag, Properties::TTag::FlowProblemTPFANoConvectiveMixing>(
        [] { return !deckNeedsConvectiveMixing(); });
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

} // namespace Opm
