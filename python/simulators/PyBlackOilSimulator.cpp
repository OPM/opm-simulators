/*
  Copyright 2020 Equinor ASA.

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
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <flow/flow_ebos_blackoil.hpp>
// NOTE: EXIT_SUCCESS, EXIT_FAILURE is defined in cstdlib
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <opm/simulators/flow/python/PyBlackOilSimulator.hpp>

namespace py = pybind11;

namespace Opm::Pybind {
PyBlackOilSimulator::PyBlackOilSimulator( const std::string &deckFilename)
    : deckFilename_{deckFilename}
{
}

PyBlackOilSimulator::PyBlackOilSimulator(
    std::shared_ptr<Opm::Deck> deck,
    std::shared_ptr<Opm::EclipseState> state,
    std::shared_ptr<Opm::Schedule> schedule,
    std::shared_ptr<Opm::SummaryConfig> summary_config
)
    : deck_{std::move(deck)}
    , eclipse_state_{std::move(state)}
    , schedule_{std::move(schedule)}
    , summary_config_{std::move(summary_config)}
{
}

bool PyBlackOilSimulator::checkSimulationFinished()
{
    return this->mainEbos_->getSimTimer()->done();
}

const Opm::FlowMainEbos<typename Opm::Pybind::PyBlackOilSimulator::TypeTag>&
         PyBlackOilSimulator::getFlowMainEbos() const
{
    if (this->mainEbos_) {
        return *this->mainEbos_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMainEbos object" );
    }
}

py::array_t<double> PyBlackOilSimulator::getPorosity()
{
    std::size_t len;
    auto array = materialState_->getPorosity(&len);
    return py::array(len, array.get());
}

int PyBlackOilSimulator::run()
{
    auto mainObject = Opm::Main( deckFilename_ );
    return mainObject.runStatic<Opm::Properties::TTag::EclFlowProblemTPFA>();
}

void PyBlackOilSimulator::setPorosity( py::array_t<double,
    py::array::c_style | py::array::forcecast> array)
{
    std::size_t size_ = array.size();
    const double *poro = array.data();
    materialState_->setPorosity(poro, size_);
}

void PyBlackOilSimulator::advance(int report_step)
{
    while (currentStep() < report_step) {
        step();
    }
}

int PyBlackOilSimulator::step()
{
    if (!hasRunInit_) {
        throw std::logic_error("step() called before step_init()");
    }
    if (hasRunCleanup_) {
        throw std::logic_error("step() called after step_cleanup()");
    }
    if(checkSimulationFinished()) {
        throw std::logic_error("step() called, but simulation is done");
    }
    //if (this->debug_)
    //    this->mainEbos_->getSimTimer()->report(std::cout);
    auto result = mainEbos_->executeStep();
    return result;
}

// This returns the report step number that will be executed next time step()
//   is called.
int PyBlackOilSimulator::currentStep()
{
    return this->mainEbos_->getSimTimer()->currentStepNum();
    // NOTE: this->ebosSimulator_->episodeIndex() would also return the current
    // report step number, but this number is always delayed by 1 step relative
    // to this->mainEbos_->getSimTimer()->currentStepNum()
    // See details in runStep() in file SimulatorFullyImplicitBlackoilEbos.hpp
}


int PyBlackOilSimulator::stepCleanup()
{
    hasRunCleanup_ = true;
    return mainEbos_->executeStepsCleanup();
}

int PyBlackOilSimulator::stepInit()
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
    if (this->deck_) {
        main_ = std::make_unique<Opm::Main>(
            this->deck_->getDataFile(),
            this->eclipse_state_,
            this->schedule_,
            this->summary_config_
        );
    }
    else {
        main_ = std::make_unique<Opm::Main>( deckFilename_ );
    }
    int exitCode = EXIT_SUCCESS;
    mainEbos_ = main_->initFlowEbosBlackoil(exitCode);
    if (mainEbos_) {
        int result = mainEbos_->executeInitStep();
        hasRunInit_ = true;
        ebosSimulator_ = mainEbos_->getSimulatorPtr();
        materialState_ = std::make_unique<PyMaterialState<TypeTag>>(
            ebosSimulator_);
        return result;
    }
    else {
        return exitCode;
    }
}

void export_PyBlackOilSimulator(py::module& m)
{
    py::class_<PyBlackOilSimulator>(m, "BlackOilSimulator")
        .def(py::init< const std::string& >())
        .def(py::init<
            std::shared_ptr<Opm::Deck>,
            std::shared_ptr<Opm::EclipseState>,
            std::shared_ptr<Opm::Schedule>,
            std::shared_ptr<Opm::SummaryConfig> >())
        .def("get_porosity", &PyBlackOilSimulator::getPorosity,
            py::return_value_policy::copy)
        .def("run", &PyBlackOilSimulator::run)
        .def("set_porosity", &PyBlackOilSimulator::setPorosity)
        .def("current_step", &PyBlackOilSimulator::currentStep)
        .def("step", &PyBlackOilSimulator::step)
        .def("advance", &PyBlackOilSimulator::advance, py::arg("report_step"))
        .def("step_init", &PyBlackOilSimulator::stepInit)
        .def("step_cleanup", &PyBlackOilSimulator::stepCleanup);
}

} // namespace Opm::Pybind

