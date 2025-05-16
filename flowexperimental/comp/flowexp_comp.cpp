/*
  Copyright 2024, SINTEF Digital

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

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <fmt/format.h>

#include <array>
#include <cstdlib>
#include <tuple>

#include "flowexp_comp.hpp"

template <int compileTimeComponent>
std::tuple<bool, int>
runComponent(int runtimeComponent, bool water, int argc, char** argv)
{
    if (runtimeComponent == compileTimeComponent) {
        if (water)
            return std::make_tuple(true, Opm::dispatchFlowExpComp<compileTimeComponent, true>(argc, argv));
        else
            return std::make_tuple(true, Opm::dispatchFlowExpComp<compileTimeComponent, false>(argc, argv));
    }
    return std::make_tuple(false, EXIT_FAILURE);
}


/**
 * @brief Runs a specified runtime component.
 *
 * This function checks if the provided runtime component matches the current compile-time component.
 * If they match, it dispatches the simulator for the number of components and returns
 * a tuple where the second element is the result of the execution.
 *
 * Otherwise, it recursively calls itself with the next component in the list.
 *
 * @param runtimecomponent The runtime component identifier to be executed.
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return A tuple containing a boolean indicating if the component was executed and an integer result of the execution.
 *
 * @note We have two non-variadic templates to be able to properly overload for the base case.
 */
template <int currentCompileTimeComponent, int nextComponent, int... components>
std::tuple<bool, int>
runComponent(int runtimecomponent, bool water, int argc, char** argv)
{
    if (currentCompileTimeComponent == runtimecomponent) {
        if (water)
            return std::make_tuple(true, Opm::dispatchFlowExpComp<currentCompileTimeComponent, true>(argc, argv));
        else
            return std::make_tuple(true, Opm::dispatchFlowExpComp<currentCompileTimeComponent, false>(argc, argv));
    }
    return runComponent<nextComponent, components...>(runtimecomponent, water, argc, argv);
}

int
main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowExpCompProblem<0, true>;
    Opm::registerEclTimeSteppingParameters<double>();

    // At the moment, this is probably as optimal as can be.
    // We only read the RUNSPEC of the Deck file to get the numComp,
    // and for this we need to first read the CLI arguments.
    Opm::setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), true, false, true, 0);

    auto inputFilename
        = Opm::FlowGenericVanguard::canonicalDeckPath(Opm::Parameters::Get<Opm::Parameters::EclDeckFileName>());

    // Only read the RUNSPEC section of the deck
    const auto deck
        = Opm::Parser {}.parseFile(inputFilename, Opm::ParseContext {}, std::vector {Opm::Ecl::SectionType::RUNSPEC});
    const auto runspec = Opm::Runspec(deck);
    const auto numComps = runspec.numComps();
    const auto& phases = runspec.phases();
    const auto wat = phases.active(Opm::Phase::WATER);

    auto [componentSupported, executionStatus]
        = runComponent<OPM_COMPILE_COMPONENTS_TEMPLATE_LIST>(numComps, wat, argc, argv);

    if (!componentSupported) {
        fmt::print("Deck has {} components, not supported. In this build of the simulator, we support the "
                   "following number of components:\n\n\t{}\n\nNote that the supported components can be changed "
                   "when configuring CMake through the OPM_COMPILE_COMPONENTS option. That is, run cmake as "
                   "eg\n\n\tcmake ... -DOPM_COMPILE_COMPONENTS=\"4;7\"\n\nto compile only components 4 and 7, "
                   "or\n\n\tcmake ... -DOPM_COMPILE_COMPONENTS=3\n\nto compile only component 3, then "
                   "recompile.\n\nExiting.\n",
                   numComps,
                   fmt::join(std::array {OPM_COMPILE_COMPONENTS_TEMPLATE_LIST}, ", "));
    }
    return executionStatus;
}
