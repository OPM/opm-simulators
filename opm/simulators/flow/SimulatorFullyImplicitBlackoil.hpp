/*
  Copyright 2013, 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 Andreas Lauser
  Copyright 2017 IRIS

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

#ifndef OPM_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_HEADER_INCLUDED
#define OPM_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingEnabled.hpp>

#ifdef RESERVOIR_COUPLING_ENABLED
#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>
#include <opm/common/Exceptions.hpp>
#endif

#include <opm/grid/utility/StopWatch.hpp>

#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/flow/BlackoilModel.hpp>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>
#include <opm/simulators/flow/ExtraConvergenceOutputThread.hpp>
#include <opm/simulators/flow/NonlinearSolver.hpp>
#include <opm/simulators/flow/SimulatorConvergenceOutput.hpp>
#include <opm/simulators/flow/SimulatorReportBanners.hpp>
#include <opm/simulators/flow/SimulatorSerializer.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/wells/WellState.hpp>

#if HAVE_HDF5
#include <opm/simulators/utils/HDF5Serializer.hpp>
#endif

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace Opm::Parameters {

struct EnableAdaptiveTimeStepping { static constexpr bool value = true; };
struct OutputExtraConvergenceInfo { static constexpr auto* value = "none"; };
struct SaveStep { static constexpr auto* value = ""; };
struct SaveFile { static constexpr auto* value = ""; };
struct LoadFile { static constexpr auto* value = ""; };
struct LoadStep { static constexpr int value = -1; };
struct Slave { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm::detail {

void registerSimulatorParameters();

}

namespace Opm {

/// a simulator for the blackoil model
template<class TypeTag>
class SimulatorFullyImplicitBlackoil : private SerializableSim
{
protected:
    struct MPI_Comm_Deleter;
public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using BlackoilIndices = GetPropType<TypeTag, Properties::Indices>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using AquiferModel = GetPropType<TypeTag, Properties::AquiferModel>;
    using Model = GetPropType<TypeTag, Properties::NonlinearSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using TimeStepper = AdaptiveTimeStepping<TypeTag>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using BioeffectsModule = BlackOilBioeffectsModule<TypeTag>;

    using Solver = NonlinearSolver<TypeTag, Model>;
    using ModelParameters = typename Model::ModelParameters;
    using SolverParameters = typename Solver::SolverParameters;
    using WellModel = BlackoilWellModel<TypeTag>;

    /// Initialise from parameters and objects to observe.
    /// \param simulator Reference to main simulator
    explicit SimulatorFullyImplicitBlackoil(Simulator& simulator);

    ~SimulatorFullyImplicitBlackoil() override;

    static void registerParameters();

    /// Run the simulation.
    /// This will run succesive timesteps until timer.done() is true. It will
    /// modify the reservoir and well states.
    /// \param[in,out] timer       governs the requested reporting timesteps
    /// \param[in,out] state       state of reservoir: pressure, fluxes
    /// \return                    simulation report, with timing data
#ifdef RESERVOIR_COUPLING_ENABLED
    SimulatorReport run(SimulatorTimer& timer, int argc, char** argv);

    // This method should only be called if slave mode (i.e. Parameters::Get<Parameters::Slave>())
    // is false. We try to determine if this is a normal flow simulation or a reservoir
    // coupling master. It is a normal flow simulation if the schedule does not contain
    // any SLAVES and GRUPMAST keywords.
    bool checkRunningAsReservoirCouplingMaster();

    // NOTE: The argc and argv will be used when launching a slave process
    void init(const SimulatorTimer& timer, int argc, char** argv);
#else
    SimulatorReport run(SimulatorTimer& timer);

    void init(const SimulatorTimer& timer);
#endif

    void updateTUNING(const Tuning& tuning);

    void updateTUNINGDP(const TuningDp& tuning_dp);

    bool runStep(SimulatorTimer& timer);

    SimulatorReport finalize();

    const Grid& grid() const { return simulator_.vanguard().grid(); }

    template<class Serializer>
    void serializeOp(Serializer& serializer);

    const Model& model() const { return solver_->model(); }

protected:
    //! \brief Load simulator state from hdf5 serializer.
    void loadState(HDF5Serializer& serializer, const std::string& groupName) override;

    //! \brief Save simulator state using hdf5 serializer.
    void saveState(HDF5Serializer& serializer, const std::string& groupName) const override;

    //! \brief Returns header data
    std::array<std::string,5> getHeader() const override;

    //! \brief Returns local-to-global cell mapping.
    const std::vector<int>& getCellMapping() const override {
        return simulator_.vanguard().globalCell();
    }

    std::unique_ptr<Solver> createSolver(WellModel& wellModel);

    const EclipseState& eclState() const { return simulator_.vanguard().eclState(); }

    const Schedule& schedule() const { return simulator_.vanguard().schedule(); }

    bool isRestart() const { return eclState().getInitConfig().restartRequested(); }

    WellModel& wellModel_() { return simulator_.problem().wellModel(); }

    const WellModel& wellModel_() const { return simulator_.problem().wellModel(); }

    // Data.
    Simulator& simulator_;

    ModelParameters modelParam_;
    SolverParameters solverParam_;

    std::unique_ptr<Solver> solver_;

    // Observed objects.
    // Misc. data
    bool terminalOutput_;

    SimulatorReport report_;
    std::unique_ptr<time::StopWatch> solverTimer_;
    std::unique_ptr<time::StopWatch> totalTimer_;
    std::unique_ptr<TimeStepper> adaptiveTimeStepping_;

    SimulatorConvergenceOutput convergence_output_{};

#ifdef RESERVOIR_COUPLING_ENABLED
    bool slaveMode_{false};
    std::unique_ptr<ReservoirCouplingMaster<Scalar>> reservoirCouplingMaster_{nullptr};
    std::unique_ptr<ReservoirCouplingSlave<Scalar>> reservoirCouplingSlave_{nullptr};
#endif

    SimulatorSerializer serializer_;
};

} // namespace Opm

#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil_impl.hpp>

#endif // OPM_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_HEADER_INCLUDED
