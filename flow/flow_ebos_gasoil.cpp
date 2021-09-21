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

#include <flow/flow_ebos_gasoil.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowGasOilProblem {
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}

//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::EclFlowGasOilProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::EclFlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    typedef BlackOilTwoPhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                    getPropValue<TypeTag, Properties::EnableExtbo>(),
                                    getPropValue<TypeTag, Properties::EnablePolymer>(),
                                    getPropValue<TypeTag, Properties::EnableEnergy>(),
                                    getPropValue<TypeTag, Properties::EnableFoam>(),
                                    getPropValue<TypeTag, Properties::EnableBrine>(),
                                    /*PVOffset=*/0,
                                    /*disabledCompIdx=*/FluidSystem::waterCompIdx> type;
};
}}

namespace Opm {
void flowEbosGasOilSetDeck(double setupTime, std::shared_ptr<Deck> deck,
                           std::shared_ptr<EclipseState> eclState,
                           std::shared_ptr<Schedule> schedule,
                           std::shared_ptr<SummaryConfig> summaryConfig)
{
    using TypeTag = Properties::TTag::EclFlowGasOilProblem;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;

    Vanguard::setExternalSetupTime(setupTime);
    Vanguard::setExternalDeck(deck);
    Vanguard::setExternalEclState(eclState);
    Vanguard::setExternalSchedule(schedule);
    Vanguard::setExternalSummaryConfig(summaryConfig);
}


// ----------------- Main program -----------------
int flowEbosGasOilMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    FlowMainEbos<Properties::TTag::EclFlowGasOilProblem>
        mainfunc {argc, argv, outputCout, outputFiles} ;
    return mainfunc.execute();
}

}
