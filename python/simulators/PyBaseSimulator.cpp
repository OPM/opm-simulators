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
#include <opm/simulators/flow/python/PyBaseSimulator.hpp>
#include <opm/simulators/flow/TTagFlowProblemGasWater.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>

namespace Opm {
template<class TypeTag>
std::unique_ptr<FlowMain<TypeTag>>
flowMainInit(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    return std::make_unique<FlowMain<TypeTag>>(argc, argv, outputCout, outputFiles);
}

// Explicit instantiation of the above flowMainInit() for the two problem types
template std::unique_ptr<FlowMain<Properties::TTag::FlowProblemTPFA>>
flowMainInit<Properties::TTag::FlowProblemTPFA>(
    int argc, char** argv, bool outputCout, bool outputFiles
);
template std::unique_ptr<FlowMain<Properties::TTag::FlowGasWaterProblem>>
flowMainInit<Properties::TTag::FlowGasWaterProblem>(
    int argc, char** argv, bool outputCout, bool outputFiles
);

}  // namespace Opm

namespace py = pybind11;

namespace Opm::Pybind {
template<class TypeTag>
PyBaseSimulator<TypeTag>::PyBaseSimulator(const std::string& deck_filename,
                                          const std::vector<std::string>& args)
    : deck_filename_{deck_filename}
    , args_{args}
{
}

template<class TypeTag>
PyBaseSimulator<TypeTag>::PyBaseSimulator(
    std::shared_ptr<Opm::Deck> deck,
    std::shared_ptr<Opm::EclipseState> state,
    std::shared_ptr<Opm::Schedule> schedule,
    std::shared_ptr<Opm::SummaryConfig> summary_config,
    const std::vector<std::string>& args
)
    : deck_{std::move(deck)}
    , eclipse_state_{std::move(state)}
    , schedule_{std::move(schedule)}
    , summary_config_{std::move(summary_config)}
    , args_{args}
{
}

// Public methods alphabetically sorted
// ------------------------------------
template<class TypeTag>
void PyBaseSimulator<TypeTag>::advance(int report_step)
{
    while (currentStep() < report_step) {
        step();
    }
}

template<class TypeTag>
bool PyBaseSimulator<TypeTag>::checkSimulationFinished()
{
    return getFlowMain().getSimTimer()->done();
}

// This returns the report step number that will be executed next time step()
//   is called.
template<class TypeTag>
int PyBaseSimulator<TypeTag>::currentStep()
{
    return getFlowMain().getSimTimer()->currentStepNum();
    // NOTE: this->simulator_->episodeIndex() would also return the current
    // report step number, but this number is always delayed by 1 step relative
    // to this->flow_main_->getSimTimer()->currentStepNum()
    // See details in runStep() in file SimulatorFullyImplicitBlackoilEbos.hpp
}

template<class TypeTag>
py::array_t<double> PyBaseSimulator<TypeTag>::getCellVolumes() {
    auto vector = getMaterialState().getCellVolumes();
    return py::array(vector.size(), vector.data());
}

template<class TypeTag>
double PyBaseSimulator<TypeTag>::getDT() {
    return getFlowMain().getPreviousReportStepSize();
}

template<class TypeTag>
py::array_t<double> PyBaseSimulator<TypeTag>::getPorosity()
{
    auto vector = getMaterialState().getPorosity();
    return py::array(vector.size(), vector.data());
}

template<class TypeTag>
py::array_t<double>
PyBaseSimulator<TypeTag>::
getFluidStateVariable(const std::string &name) const
{
    auto vector = getFluidState().getFluidStateVariable(name);
    return py::array(vector.size(), vector.data());
}

template<class TypeTag>
py::array_t<double>
PyBaseSimulator<TypeTag>::
getPrimaryVariable(const std::string &variable) const
{
    auto vector = getFluidState().getPrimaryVariable(variable);
    return py::array(vector.size(), vector.data());
}

template<class TypeTag>
py::array_t<int>
PyBaseSimulator<TypeTag>::
getPrimaryVarMeaning(const std::string &variable) const
{
    auto vector = getFluidState().getPrimaryVarMeaning(variable);
    return py::array(vector.size(), vector.data());
}

template<class TypeTag>
std::map<std::string, int>
PyBaseSimulator<TypeTag>::
getPrimaryVarMeaningMap(const std::string &variable) const
{

    return getFluidState().getPrimaryVarMeaningMap(variable);
}

// template<class TypeTag>
// int PyBaseSimulator<TypeTag>::run()
// {
//     auto main_object = Opm::Main( this->deck_filename_ );
//     return main_object.runStatic<Opm::Properties::TTag::FlowGasWaterProblem>();
// }

template<class TypeTag>
void PyBaseSimulator<TypeTag>::setPorosity( py::array_t<double,
    py::array::c_style | py::array::forcecast> array)
{
    std::size_t size_ = array.size();
    const double *poro = array.data();
    getMaterialState().setPorosity(poro, size_);
}

template<class TypeTag>
void
PyBaseSimulator<TypeTag>::
setPrimaryVariable(
    const std::string &variable,
    py::array_t<double,
    py::array::c_style | py::array::forcecast> array
)
{
    std::size_t size_ = array.size();
    const double *data = array.data();
    getFluidState().setPrimaryVariable(variable, data, size_);
}

template<class TypeTag>
void PyBaseSimulator<TypeTag>::setupMpi(bool mpi_init, bool mpi_finalize)
{
    if (this->has_run_init_) {
        throw std::logic_error("mpi_init() called after step_init()");
    }
    this->mpi_init_ = mpi_init;
    this->mpi_finalize_ = mpi_finalize;
}

template<class TypeTag>
int PyBaseSimulator<TypeTag>::step()
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

template<class TypeTag>
int PyBaseSimulator<TypeTag>::stepCleanup()
{
    this->has_run_cleanup_ = true;
    return getFlowMain().executeStepsCleanup();
}

template<class TypeTag>
int PyBaseSimulator<TypeTag>::stepInit()
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
        this->main_ = std::make_unique<Opm::PyMain<TypeTag>>(
            this->deck_->getDataFile(),
            this->eclipse_state_,
            this->schedule_,
            this->summary_config_,
            this->mpi_init_,
            this->mpi_finalize_
        );
    }
    else {
        this->main_ = std::make_unique<Opm::PyMain<TypeTag>>(
            this->deck_filename_,
            this->mpi_init_,
            this->mpi_finalize_
        );
    }
    this->main_->setArguments(args_);
    int exit_code = EXIT_SUCCESS;
    this->flow_main_ = this->main_->initFlowBlackoil(exit_code);
    if (this->flow_main_) {
        int result = this->flow_main_->executeInitStep();
        this->has_run_init_ = true;
        this->simulator_ = this->flow_main_->getSimulatorPtr();
        this->fluid_state_ = std::make_unique<PyFluidState<TypeTag>>(this->simulator_);
        this->material_state_ = std::make_unique<PyMaterialState<TypeTag>>(this->simulator_);
        return result;
    }
    else {
        return exit_code;
    }
}

template<class TypeTag>
int PyBaseSimulator<TypeTag>::run()
{
    auto main_object = Opm::Main( this->deck_filename_ );
    return main_object.runStatic<TypeTag>();
}

// Private methods
// ---------------
template<class TypeTag>
Opm::FlowMain<TypeTag>&
PyBaseSimulator<TypeTag>::getFlowMain() const
{
    if (this->flow_main_) {
        return *this->flow_main_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMain object" );
    }
}

template<class TypeTag>
PyFluidState<TypeTag>&
PyBaseSimulator<TypeTag>::getFluidState() const
{
    if (this->fluid_state_) {
        return *this->fluid_state_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMainEbos object" );
    }
}

template<class TypeTag>
PyMaterialState<TypeTag>&
PyBaseSimulator<TypeTag>::getMaterialState() const
{
    if (this->material_state_) {
        return *this->material_state_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMain object" );
    }
}

// NOTE: We need the below explicit instantiations or else the the symbols
//       will not be available in the shared library and we will get
//       undefined symbol errors when trying to import the module in Python.
template class PyBaseSimulator<Opm::Properties::TTag::FlowProblemTPFA>;
template class PyBaseSimulator<Opm::Properties::TTag::FlowGasWaterProblem>;

}  // namespace Opm::Pybind
