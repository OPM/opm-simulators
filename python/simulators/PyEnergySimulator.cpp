
#include "config.h"
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/simulators/flow/python/PyEnergySimulator.hpp>
#include <opm/simulators/flow/python/PyFluidState.hpp>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#include <flow/flow_ebos_energy.hpp>
// NOTE: EXIT_SUCCESS, EXIT_FAILURE is defined in cstdlib
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace Opm::Pybind {
PyEnergySimulator::
PyEnergySimulator( const std::string &deck_filename)
    : deck_filename_{deck_filename}
{
}

PyEnergySimulator::
PyEnergySimulator(
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

//used to advance the state of the simulator until the specified reporting step is reached.
void
PyEnergySimulator::
advance(int report_step)
{
    while (currentStep() < report_step) {
        step();
    }
}

bool
PyEnergySimulator::
checkSimulationFinished() const
{
    return getFlowMainEbos().getSimTimer()->done();
}

// This returns the report step number that will be executed next time step()
//   is called.
int
PyEnergySimulator::
currentStep() const
{
    return getFlowMainEbos().getSimTimer()->currentStepNum();
    // NOTE: this->ebos_simulator_->episodeIndex() would also return the current
    // report step number, but this number is always delayed by 1 step relative
    // to this->main_ebos_->getSimTimer()->currentStepNum()
    // See details in runStep() in file SimulatorFullyImplicitBlackoilEbos.hpp
}

py::array_t<double>
PyEnergySimulator::
getCellVolumes() const
{
    std::size_t len;
    auto array = getMaterialState().getCellVolumes(&len);
    return py::array(len, array.get());
}

double
PyEnergySimulator::
getDT() const
{
    return getFlowMainEbos().getPreviousReportStepSize();
}

py::array_t<double>
PyEnergySimulator::
getPorosity() const
{
    std::size_t len;
    auto array = getMaterialState().getPorosity(&len);
    return py::array(len, array.get());
}

py::array_t<double>
PyEnergySimulator::
getFluidStateVariable(const std::string &name) const
{
    std::size_t len;  //size_t unsigned int
    auto array = getFluidState().getFluidStateVariable(name, &len);
    return py::array(len, array.get());
}

py::array_t<double>
PyEnergySimulator::
getPrimaryVariable(const std::string &variable) const
{
    std::size_t len;
    auto array = getFluidState().getPrimaryVariable(variable, &len);
    return py::array(len, array.get());
}

py::array_t<int>
PyEnergySimulator::
getPrimaryVarMeaning(const std::string &variable) const
{
    std::size_t len;
    auto array = getFluidState().getPrimaryVarMeaning(variable, &len);
    return py::array(len, array.get());
}

std::map<std::string, int>
PyEnergySimulator::
getPrimaryVarMeaningMap(const std::string &variable) const
{

    return getFluidState().getPrimaryVarMeaningMap(variable);
}

int
PyEnergySimulator::
run()
{
    auto main_object = ::Opm::Main( this->deck_filename_ );
    return main_object.runStatic<::Opm::Properties::TTag::EclEnergyProblemTPFA>();
}

void
PyEnergySimulator::
setPorosity( py::array_t<double, py::array::c_style | py::array::forcecast> array)
{
    std::size_t size_ = array.size();
    const double *poro = array.data();
    getMaterialState().setPorosity(poro, size_);
}

void
PyEnergySimulator::
setPrimaryVariable(
    const std::string &idx_name,
    py::array_t<double,py::array::c_style | py::array::forcecast> array
)
{
    std::size_t size_ = array.size();
    const double *data = array.data();
    getFluidState().setPrimaryVariable(idx_name, data, size_);
}

int
PyEnergySimulator::
step()
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
    auto result = getFlowMainEbos().executeStep();
    return result;
}

int
PyEnergySimulator::
stepCleanup()
{
    this->has_run_cleanup_ = true;
    return getFlowMainEbos().executeStepsCleanup();
}

int
PyEnergySimulator::
stepInit()
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
        this->main_ = std::make_unique<::Opm::Main>(
            this->deck_->getDataFile(),
            this->eclipse_state_,
            this->schedule_,
            this->summary_config_
        );
    }
    else {
        this->main_ = std::make_unique<::Opm::Main>( this->deck_filename_ );
    }
    int exit_code = EXIT_SUCCESS;
    this->main_ebos_ = this->main_->initFlowEbosEnergy(exit_code);
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

::Opm::FlowMainEbos<typename PyEnergySimulator::TypeTag>&
PyEnergySimulator::
getFlowMainEbos() const
{
    if (this->main_ebos_) {
        return *this->main_ebos_;
    }
    else {
        throw std::runtime_error("FlowSimulator not initialized: "
            "Cannot get reference to FlowMainEbos object" );
    }
}

PyFluidState<typename PyEnergySimulator::TypeTag>&
PyEnergySimulator::
getFluidState() const
{
    if (this->fluid_state_) {
        return *this->fluid_state_;
    }
    else {
        throw std::runtime_error("FlowSimulator not initialized: "
            "Cannot get reference to FlowMainEbos object" );
    }
}

PyMaterialState<typename Opm::Pybind::PyEnergySimulator::TypeTag>&
PyEnergySimulator::
getMaterialState() const
{
    if (this->material_state_) {
        return *this->material_state_;
    }
    else {
        throw std::runtime_error("FlowSimulator not initialized: "
            "Cannot get reference to FlowMainEbos object" );
    }
}

// Exported functions
void export_PyEnergySimulator(py::module& m)
{
    py::class_<PyEnergySimulator>(m, "EnergySimulator")
        .def(py::init< const std::string& >())
        .def(py::init<
            std::shared_ptr<Opm::Deck>,
            std::shared_ptr<Opm::EclipseState>,
            std::shared_ptr<Opm::Schedule>,
            std::shared_ptr<Opm::SummaryConfig> >())
        .def("advance", &PyEnergySimulator::advance, py::arg("report_step"))
        .def("current_step", &PyEnergySimulator::currentStep)
        .def("get_cell_volumes", &PyEnergySimulator::getCellVolumes,
            py::return_value_policy::copy)
        .def("get_dt", &PyEnergySimulator::getDT)
        .def("get_fluidstate_variable", &PyEnergySimulator::getFluidStateVariable,
            py::return_value_policy::copy, py::arg("name"))
        .def("get_porosity", &PyEnergySimulator::getPorosity,
            py::return_value_policy::copy)
        .def("get_primary_variable_meaning", &PyEnergySimulator::getPrimaryVarMeaning,
            py::return_value_policy::copy, py::arg("variable"))
        .def("get_primary_variable_meaning_map", &PyEnergySimulator::getPrimaryVarMeaningMap,
            py::return_value_policy::copy, py::arg("variable"))
        .def("get_primary_variable", &PyEnergySimulator::getPrimaryVariable,
            py::return_value_policy::copy, py::arg("variable"))
        .def("run", &PyEnergySimulator::run)
        .def("set_porosity", &PyEnergySimulator::setPorosity)
        .def("set_primary_variable", &PyEnergySimulator::setPrimaryVariable,
            py::arg("idx_name"), py::arg("value"))
        .def("step", &PyEnergySimulator::step)
        .def("step_cleanup", &PyEnergySimulator::stepCleanup)
        .def("step_init", &PyEnergySimulator::stepInit);
}

} // namespace Opm::Pybind

