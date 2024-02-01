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
#include <opm/simulators/flow/FlowMain.hpp>
//#include <opm/simulators/flow/python/PyFluidState.hpp>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#include <opm/simulators/flow/python/PyBlackOilSimulator.hpp>
#include <opm/simulators/flow/python/PyFluidState.hpp>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#include <flow/flow_ebos_blackoil.hpp>
// NOTE: EXIT_SUCCESS, EXIT_FAILURE is defined in cstdlib
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <stdlib.h>

namespace py = pybind11;

namespace Opm::Pybind {
PyBlackOilSimulator::PyBlackOilSimulator( const std::string &deck_filename)
    : deck_filename_{deck_filename}
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

// Public methods alphabetically sorted
// ------------------------------------

void PyBlackOilSimulator::advance(int report_step)
{
    while (currentStep() < report_step) {
        step();
    }
}

bool PyBlackOilSimulator::checkSimulationFinished()
{
    return getFlowMain().getSimTimer()->done();
}

// This returns the report step number that will be executed next time step()
//   is called.
int PyBlackOilSimulator::currentStep()
{
    return getFlowMain().getSimTimer()->currentStepNum();
    // NOTE: this->ebos_simulator_->episodeIndex() would also return the current
    // report step number, but this number is always delayed by 1 step relative
    // to this->main_ebos_->getSimTimer()->currentStepNum()
    // See details in runStep() in file SimulatorFullyImplicitBlackoilEbos.hpp
}

py::array_t<double> PyBlackOilSimulator::getCellVolumes() {
    auto vector = getMaterialState().getCellVolumes();
    return py::array(vector.size(), vector.data());
}

double PyBlackOilSimulator::getDT() {
    return getFlowMain().getPreviousReportStepSize();
}

py::array_t<double> PyBlackOilSimulator::getPorosity()
{
    auto vector = getMaterialState().getPorosity();
    return py::array(vector.size(), vector.data());
}

py::array_t<double>
PyBlackOilSimulator::
getFluidStateVariable(const std::string &name) const
{
    setenv("OPM_DEBUG", "1", /*overwrite=*/1);
    std::cout << "getFluidStateVariable: " << name << std::endl;
    auto vector = getFluidState().getFluidStateVariable(name);
    return py::array(vector.size(), vector.data());
}

py::array_t<double>
PyBlackOilSimulator::
getPrimaryVariable(const std::string &variable) const
{
    auto vector = getFluidState().getPrimaryVariable(variable);
    return py::array(vector.size(), vector.data());
}

py::array_t<int>
PyBlackOilSimulator::
getPrimaryVarMeaning(const std::string &variable) const
{
    auto vector = getFluidState().getPrimaryVarMeaning(variable);
    return py::array(vector.size(), vector.data());
}

std::map<std::string, int>
PyBlackOilSimulator::
getPrimaryVarMeaningMap(const std::string &variable) const
{
    return getFluidState().getPrimaryVarMeaningMap(variable);
}

int PyBlackOilSimulator::run()
{
    auto main_object = Opm::Main( this->deck_filename_ );
    return main_object.runStatic<Opm::Properties::TTag::EclFlowProblemTPFA>();
}

void PyBlackOilSimulator::setPorosity( py::array_t<double,
    py::array::c_style | py::array::forcecast> array)
{
    std::size_t size_ = array.size();
    const double *poro = array.data();
    getMaterialState().setPorosity(poro, size_);
}

void
PyBlackOilSimulator::
setPrimaryVariable(
    const std::string &idx_name,
    py::array_t<double,
    py::array::c_style | py::array::forcecast> array
)
{
    std::size_t size_ = array.size();
    const double *data = array.data();
    getFluidState().setPrimaryVariable(idx_name, data, size_);
}

int PyBlackOilSimulator::step()
{
    if (!this->has_run_init_) {
        throw std::logic_error("step() called before step_init()");
    }
    if (this->has_run_cleanup_) {
        throw std::logic_error("step() called after step_cleanup()");
    }
    if(checkSimulationFinished()) {
        throw std::logic_error("step() called, but simulation is done");
    }
    //if (this->debug_)
    //    this->mainEbos_->getSimTimer()->report(std::cout);
    auto result = getFlowMain().executeStep();
    return result;
}

int PyBlackOilSimulator::stepCleanup()
{
    this->has_run_cleanup_ = true;
    return getFlowMain().executeStepsCleanup();
}

int PyBlackOilSimulator::stepInit()
{

    if (this->has_run_init_) {
        // Running step_init() multiple times is not implemented yet,
        if (this->has_run_cleanup_) {
            throw std::logic_error("step_init() called again");
        }
        else {
            return EXIT_SUCCESS;
        }
    }
    if (this->deck_) {
        this->main_ = std::make_unique<Opm::Main>(
            this->deck_->getDataFile(),
            this->eclipse_state_,
            this->schedule_,
            this->summary_config_
        );
    }
    else {
        this->main_ = std::make_unique<Opm::Main>( this->deck_filename_ );
    }
    int exit_code = EXIT_SUCCESS;
    this->main_ebos_ = this->main_->initFlowEbosBlackoil(exit_code);
    if (this->main_ebos_) {
        int result = this->main_ebos_->executeInitStep();
        this->has_run_init_ = true;
        this->ebos_simulator_ = this->main_ebos_->getSimulatorPtr();
        this->fluid_state_ = std::make_unique<PyFluidState<TypeTag>>(this->ebos_simulator_);
        this->material_state_ = std::make_unique<PyMaterialState<TypeTag>>(this->ebos_simulator_);
        return result;
    }
    else {
        return exit_code;
    }
}

// Private methods alphabetically sorted
// ------------------------------------

Opm::FlowMain<typename Opm::Pybind::PyBlackOilSimulator::TypeTag>&
         PyBlackOilSimulator::getFlowMain() const
{
    if (this->main_ebos_) {
        return *this->main_ebos_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMain object" );
    }
}

PyFluidState<typename PyBlackOilSimulator::TypeTag>&
PyBlackOilSimulator::
getFluidState() const
{
    if (this->fluid_state_) {
        return *this->fluid_state_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMain object" );
    }
}

PyMaterialState<typename PyBlackOilSimulator::TypeTag>&
PyBlackOilSimulator::getMaterialState() const
{
    if (this->material_state_) {
        return *this->material_state_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMain object" );
    }
}

// Exported functions
void export_PyBlackOilSimulator(py::module& m)
{
    py::class_<PyBlackOilSimulator>(m, "BlackOilSimulator")
        .def(py::init< const std::string& >())
        .def(py::init<
            std::shared_ptr<Opm::Deck>,
            std::shared_ptr<Opm::EclipseState>,
            std::shared_ptr<Opm::Schedule>,
            std::shared_ptr<Opm::SummaryConfig> >())
        .def("advance", &PyBlackOilSimulator::advance, py::arg("report_step"))
        .def("current_step", &PyBlackOilSimulator::currentStep)
        .def("get_cell_volumes", &PyBlackOilSimulator::getCellVolumes,
            py::return_value_policy::copy)
        .def("get_dt", &PyBlackOilSimulator::getDT)
        .def("get_fluidstate_variable", &PyBlackOilSimulator::getFluidStateVariable,
            py::return_value_policy::copy, py::arg("name"))
        .def("get_porosity", &PyBlackOilSimulator::getPorosity,
            py::return_value_policy::copy)
        .def("get_primary_variable_meaning", &PyBlackOilSimulator::getPrimaryVarMeaning,
            py::return_value_policy::copy, py::arg("variable"))
        .def("get_primary_variable_meaning_map", &PyBlackOilSimulator::getPrimaryVarMeaningMap,
            py::return_value_policy::copy, py::arg("variable"))
        .def("get_primary_variable", &PyBlackOilSimulator::getPrimaryVariable,
            py::return_value_policy::copy, py::arg("variable"))
        .def("run", &PyBlackOilSimulator::run)
        .def("set_porosity", &PyBlackOilSimulator::setPorosity)
        .def("set_primary_variable", &PyBlackOilSimulator::setPrimaryVariable,
            py::arg("idx_name"), py::arg("value"))
        .def("current_step", &PyBlackOilSimulator::currentStep)
        .def("step", &PyBlackOilSimulator::step)
        .def("step_cleanup", &PyBlackOilSimulator::stepCleanup)
        .def("step_init", &PyBlackOilSimulator::stepInit);
}

} // namespace Opm::Pybind

