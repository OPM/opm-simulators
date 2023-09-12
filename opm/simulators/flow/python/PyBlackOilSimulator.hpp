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

#ifndef OPM_PY_BLACKOIL_SIMULATOR_HEADER_INCLUDED
#define OPM_PY_BLACKOIL_SIMULATOR_HEADER_INCLUDED

#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/simulators/flow/python/Pybind11Exporter.hpp>
#include <opm/simulators/flow/python/PyFluidState.hpp>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <map>

namespace Opm::Pybind {
class PyBlackOilSimulator
{
private:
    using TypeTag = ::Opm::Properties::TTag::EclFlowProblemTPFA;
    using Simulator = ::Opm::GetPropType<TypeTag, ::Opm::Properties::Simulator>;

public:
    PyBlackOilSimulator( const std::string& deckFilename);
    PyBlackOilSimulator(
        std::shared_ptr<::Opm::Deck> deck,
        std::shared_ptr<::Opm::EclipseState> state,
        std::shared_ptr<::Opm::Schedule> schedule,
        std::shared_ptr<::Opm::SummaryConfig> summary_config);
    void advance(int report_step);
    bool checkSimulationFinished() const;
    int currentStep() const;
    py::array_t<double> getCellVolumes() const;
    double getDT() const;
    py::array_t<double> getFluidStateVariable(const std::string &name) const;
    py::array_t<double> getPorosity() const;
    py::array_t<double> getPrimaryVariable(const std::string &variable) const;
    py::array_t<int> getPrimaryVarMeaning(const std::string &variable) const;
    std::map<std::string, int> getPrimaryVarMeaningMap(const std::string &variable) const;
    int run();
    void setPorosity(py::array_t<double, py::array::c_style | py::array::forcecast> array);
    void setPrimaryVariable(
        const std::string &idx_name,
        py::array_t<double,
        py::array::c_style | py::array::forcecast> array);
    int step();
    int stepCleanup();
    int stepInit();

private:
    ::Opm::FlowMainEbos<TypeTag>& getFlowMainEbos() const;
    PyFluidState<TypeTag>& getFluidState() const;
    PyMaterialState<TypeTag>& getMaterialState() const;

    const std::string deck_filename_;
    bool has_run_init_ = false;
    bool has_run_cleanup_ = false;
    //bool debug_ = false;
    // This *must* be declared before other pointers
    // to simulator objects. This in order to deinitialize
    // MPI at the correct time (ie after the other objects).
    std::unique_ptr<::Opm::Main> main_;

    std::unique_ptr<::Opm::FlowMainEbos<TypeTag>> main_ebos_;
    Simulator *ebos_simulator_;
    std::unique_ptr<PyFluidState<TypeTag>> fluid_state_;
    std::unique_ptr<PyMaterialState<TypeTag>> material_state_;
    std::shared_ptr<::Opm::Deck> deck_;
    std::shared_ptr<::Opm::EclipseState> eclipse_state_;
    std::shared_ptr<::Opm::Schedule> schedule_;
    std::shared_ptr<::Opm::SummaryConfig> summary_config_;
};

} // namespace Opm::Pybind
#endif // OPM_PY_BLACKOIL_SIMULATOR_HEADER_INCLUDED
