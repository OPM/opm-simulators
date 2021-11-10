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
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
// NOTE: EXIT_SUCCESS, EXIT_FAILURE is defined in cstdlib
#include <cstdlib>
#include <iostream>
#include <string>
#include <opm/simulators/flow/python/simulators.hpp>

namespace py = pybind11;

namespace Opm::Pybind {
BlackOilSimulator::BlackOilSimulator( const std::string &deckFilename)
    : deckFilename_{deckFilename}
{
}

py::array_t<double> BlackOilSimulator::getPorosity()
{
    std::size_t len;
    auto array = materialState_->getPorosity(&len);
    return py::array(len, array.get());
}

int BlackOilSimulator::run()
{
    auto mainObject = Opm::Main( deckFilename_ );
    return mainObject.runDynamic();
}

void BlackOilSimulator::setPorosity( py::array_t<double,
    py::array::c_style | py::array::forcecast> array)
{
    std::size_t size_ = array.size();
    const double *poro = array.data();
    materialState_->setPorosity(poro, size_);
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
        ebosSimulator_ = mainEbos_->getSimulatorPtr();
        materialState_ = std::make_unique<PyMaterialState<TypeTag>>(
            ebosSimulator_);
        return result;
    }
    else {
        return exitCode;
    }
}

} // namespace Opm::Pybind

PYBIND11_MODULE(simulators, m)
{
    using namespace Opm::Pybind;
    py::class_<BlackOilSimulator>(m, "BlackOilSimulator")
        .def(py::init< const std::string& >())
        .def("get_porosity", &BlackOilSimulator::getPorosity,
            py::return_value_policy::copy)
        .def("run", &BlackOilSimulator::run)
        .def("set_porosity", &BlackOilSimulator::setPorosity)
        .def("step", &BlackOilSimulator::step)
        .def("step_init", &BlackOilSimulator::stepInit)
        .def("step_cleanup", &BlackOilSimulator::stepCleanup);
}
