#include "config.h"
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#define FLOW_BLACKOIL_ONLY
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
// NOTE: EXIT_SUCCESS, EXIT_FAILURE is defined in cstdlib
#include <cstdlib>
#include <iostream>
#include <string>
#include <opm/simulators/flow/python/simulators.hpp>

namespace py = pybind11;

namespace Opm::Pybind {
BlackOilSimulator::BlackOilSimulator( const std::string &deckFilename)
    : deckFilename_{deckFilename}, hasRunInit_{false}, hasRunCleanup_{false}
{
}

int BlackOilSimulator::run()
{
    auto mainObject = Opm::Main( deckFilename_ );
    return mainObject.runDynamic();
}

int BlackOilSimulator::step()
{
    if (!hasRunInit_) {
        throw std::logic_error("step() called before step_init()");
    }
    if (hasRunCleanup_) {
        throw std::logic_error("step() called after step_cleanup()");
    }
    return mainEbos_->executeStep();
}

int BlackOilSimulator::stepCleanup()
{
    hasRunCleanup_ = true;
    return mainEbos_->executeStepsCleanup();
}

int BlackOilSimulator::stepInit()
{

    if (hasRunInit_) {
        // Running step_init() multiple times is not implemented yet,
        if (hasRunCleanup_) {
            throw std::logic_error("step_init() called again");
        }
        else {
            return EXIT_SUCCESS;
        }
    }
    main_ = std::make_unique<Opm::Main>( deckFilename_ );
    int exitCode = EXIT_SUCCESS;
    mainEbos_ = main_->initFlowEbosBlackoil(exitCode);
    if (mainEbos_) {
        int result = mainEbos_->executeInitStep();
        hasRunInit_ = true;
        return result;
    }
    else {
        return exitCode;
    }
}

} // namespace Opm::Python

PYBIND11_MODULE(simulators, m)
{
    py::class_<Opm::Pybind::BlackOilSimulator>(m, "BlackOilSimulator")
        .def(py::init< const std::string& >())
        .def("run", &Opm::Pybind::BlackOilSimulator::run)
        .def("step", &Opm::Pybind::BlackOilSimulator::step)
        .def("step_init", &Opm::Pybind::BlackOilSimulator::stepInit)
        .def("step_cleanup", &Opm::Pybind::BlackOilSimulator::stepCleanup);
}
