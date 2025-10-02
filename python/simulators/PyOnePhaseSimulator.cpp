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
#include "config.h"
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/simulators/flow/python/PyOnePhaseSimulator.hpp>
// NOTE: This file will be generated at compile time and placed in the build directory
// See python/generate_docstring_hpp.py, and python/simulators/CMakeLists.txt for details
#include <PyOnePhaseSimulatorDoc.hpp>
// NOTE: EXIT_SUCCESS, EXIT_FAILURE is defined in cstdlib
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace Opm::Properties {

    //! The indices required by the model
    template<class TypeTag>
    struct Indices<TypeTag, TTag::FlowOnePhaseProblem>
    {
    private:
        // it is unfortunately not possible to simply use 'TypeTag' here because this leads
        // to cyclic definitions of some properties. if this happens the compiler error
        // messages unfortunately are *really* confusing and not really helpful.
        using BaseTypeTag = TTag::FlowProblem;
        using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

    public:
        using type = BlackOilOnePhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                             getPropValue<TypeTag, Properties::EnableExtbo>(),
                                             getPropValue<TypeTag, Properties::EnablePolymer>(),
                                             getPropValue<TypeTag, Properties::EnableEnergy>(),
                                             getPropValue<TypeTag, Properties::EnableFoam>(),
                                             getPropValue<TypeTag, Properties::EnableBrine>(),
                                             /*PVOffset=*/0,
                                             /*enabledCompIdx=*/FluidSystem::waterCompIdx,
                                             getPropValue<TypeTag, Properties::EnableBioeffects>()>;
    };  // struct Indices

}  // namespace Opm::Properties

// NOTE: We need the below explicit instantiations or else the symbols
//       will not be available in the shared library and we will get
//       undefined symbol errors when trying to import the module in Python.
namespace Opm::Pybind {

template class PyBaseSimulator<Opm::Properties::TTag::FlowOnePhaseProblem>;

} // namespace Opm::Pybind

// Define main function
namespace Opm {

template class PyMain<Opm::Properties::TTag::FlowOnePhaseProblem>;
template std::unique_ptr<FlowMain<Opm::Properties::TTag::FlowOnePhaseProblem>>
  flowMainInit<Opm::Properties::TTag::FlowOnePhaseProblem>(
    int argc, char** argv, bool outputCout, bool outputFiles);

}  // namespace Opm

namespace py = pybind11;

namespace Opm::Pybind {

// Exported functions
void export_PyOnePhaseSimulator(py::module& m)
{
    using namespace Opm::Pybind::DocStrings;
    using TypeTag = Opm::Properties::TTag::FlowOnePhaseProblem;

    py::class_<PyBaseSimulator<TypeTag>>(
        m,
        "_BaseSimulatorOP",
        py::module_local()
    );
    py::class_<PyOnePhaseSimulator, PyBaseSimulator<TypeTag> >(m, "OnePhaseSimulator")
        .def(py::init<const std::string&,
                      const std::vector<std::string>&>(),
             PyOnePhaseSimulator_filename_constructor_docstring,
             py::arg("filename"), py::arg("args") = std::vector<std::string>{})
        .def(py::init<
             std::shared_ptr<Opm::Deck>,
             std::shared_ptr<Opm::EclipseState>,
             std::shared_ptr<Opm::Schedule>,
             std::shared_ptr<Opm::SummaryConfig>,
             const std::vector<std::string>&>(),
             PyOnePhaseSimulator_objects_constructor_docstring,
             py::arg("Deck"), py::arg("EclipseState"), py::arg("Schedule"), py::arg("SummaryConfig"),
             py::arg("args") = std::vector<std::string>{})
        .def("advance", &PyBaseSimulator<TypeTag>::advance, advance_docstring, py::arg("report_step"))
        .def("check_simulation_finished", &PyBaseSimulator<TypeTag>::checkSimulationFinished,
             checkSimulationFinished_docstring)
        .def("current_step", &PyBaseSimulator<TypeTag>::currentStep, currentStep_docstring)
        .def("get_cell_volumes", &PyBaseSimulator<TypeTag>::getCellVolumes, getCellVolumes_docstring)
        .def("get_dt", &PyBaseSimulator<TypeTag>::getDT, getDT_docstring)
        .def("get_fluidstate_variable", &PyBaseSimulator<TypeTag>::getFluidStateVariable,
            py::return_value_policy::copy, getFluidStateVariable_docstring, py::arg("name"))
        .def("get_porosity", &PyBaseSimulator<TypeTag>::getPorosity, getPorosity_docstring)
        .def("get_primary_variable_meaning", &PyBaseSimulator<TypeTag>::getPrimaryVarMeaning,
            py::return_value_policy::copy, getPrimaryVarMeaning_docstring, py::arg("variable"))
        .def("get_primary_variable_meaning_map", &PyBaseSimulator<TypeTag>::getPrimaryVarMeaningMap,
            py::return_value_policy::copy, getPrimaryVarMeaningMap_docstring, py::arg("variable"))
        .def("get_primary_variable", &PyBaseSimulator<TypeTag>::getPrimaryVariable,
            py::return_value_policy::copy, getPrimaryVariable_docstring, py::arg("variable"))
        .def("run", &PyBaseSimulator<TypeTag>::run, run_docstring)
        .def("set_porosity", &PyBaseSimulator<TypeTag>::setPorosity, setPorosity_docstring, py::arg("array"))
        .def("set_primary_variable", &PyBaseSimulator<TypeTag>::setPrimaryVariable,
            py::arg("variable"), setPrimaryVariable_docstring, py::arg("value"))
        .def("setup_mpi", &PyBaseSimulator<TypeTag>::setupMpi, setupMpi_docstring, py::arg("init"), py::arg("finalize"))
        .def("step", &PyBaseSimulator<TypeTag>::step, step_docstring)
        .def("step_cleanup", &PyBaseSimulator<TypeTag>::stepCleanup, stepCleanup_docstring)
        .def("step_init", &PyBaseSimulator<TypeTag>::stepInit, stepInit_docstring);
}

PYBIND11_MODULE(OnePhase, m)
{
    export_PyOnePhaseSimulator(m);
}

} // namespace Opm::Pybind
