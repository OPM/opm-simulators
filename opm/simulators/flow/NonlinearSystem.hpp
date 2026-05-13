#ifndef OPM_FLOW_NONLINEAR_SYSTEM_HEADER_INCLUDED
#define OPM_FLOW_NONLINEAR_SYSTEM_HEADER_INCLUDED

#include <cstddef>

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/FlowBaseProblemProperties.hpp>
#include <opm/simulators/linalg/linalgproperties.hh>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>
#include <opm/simulators/utils/ComponentName.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <span>
#include <string_view>
#include <tuple>
#include <vector>

namespace Opm {

template <class TypeTag>
class NonlinearSystem
{
public:
    // --------- Types ---------
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using ModelParameters = BlackoilModelParameters<Scalar>;
    using WellModel = GetPropType<TypeTag, Properties::WellModel>;
    using ComponentName = ::Opm::ComponentName<FluidSystem, Indices>;

    // --------- Public methods ---------
    virtual ~NonlinearSystem() = default;

    bool isParallel() const
    { return grid_.comm().size() > 1; }

    const Simulator& simulator() const
    { return simulator_; }

    Simulator& simulator()
    { return simulator_; }

    bool terminalOutputEnabled() const
    { return terminal_output_; }

    int numPhases() const
    { return Indices::numPhases; }

    const SimulatorReportSingle& failureReport() const
    { return failureReport_; }

    const std::vector<StepReport>& stepReports() const
    { return convergence_reports_; }

    const ComponentName& compNames() const
    { return compNames_; }

    const ModelParameters& param() const
    { return param_; }

    WellModel& wellModel()
    { return well_model_; }

    const WellModel& wellModel() const
    { return well_model_; }

    void beginReportStep()
    { simulator_.problem().beginEpisode(); }

    void endReportStep()
    { simulator_.problem().endEpisode(); }

    SimulatorReportSingle assembleReservoir(const SimulatorTimerInterface& timer);

    void updateTUNING(const Tuning& tuning);

    void updateTUNINGDP(const TuningDp& tuning_dp);

    void updateSolution(const GlobalEqVector& dx);

    template <class LogFailure>
    void addReservoirConvergenceMetrics(ConvergenceReport& report,
                                        const int componentIdx,
                                        const std::string_view componentName,
                                        const std::span<const Scalar> residuals,
                                        const std::span<const ConvergenceReport::ReservoirFailure::Type> types,
                                        const std::span<const Scalar> tolerances,
                                        const Scalar maxResidualAllowed,
                                        LogFailure&& logFailure) const;

protected:
    explicit NonlinearSystem(Simulator& simulator,
                             const ModelParameters& param,
                             WellModel& wellModel,
                             const bool terminal_output);

    virtual void initialLinearization(SimulatorReportSingle& report,
                                      int minIter,
                                      int maxIter,
                                      const SimulatorTimerInterface& timer);

    virtual bool shouldStoreSolutionUpdate() const
    { return false; }

    virtual void prepareSolutionUpdate()
    {}

    virtual void storeSolutionUpdate(const GlobalEqVector&)
    {}

    SimulatorReportSingle prepareStep(const SimulatorTimerInterface& timer);

    template <class WellModelType>
    SimulatorReportSingle assembleReservoir(WellModelType& wellModel);

    template <class ModelParametersType>
    void applyTUNING(ModelParametersType& param,
                     const Tuning& tuning);

    template <class ModelParametersType>
    void applyTUNINGDP(ModelParametersType& param,
                       const TuningDp& tuning_dp);

    template <class ValueType>
    std::tuple<ValueType, ValueType>
    convergenceReduction(Parallel::Communication comm,
                         const ValueType primaryVolumeLocal,
                         const ValueType secondaryVolumeLocal,
                         std::vector<ValueType>& sumValues,
                         std::vector<ValueType>& maxValues,
                         std::vector<ValueType>& averagedValues);

    void popLastStepReport()
    { convergence_reports_.back().report.pop_back(); }

    // --------- Data members ---------
    Simulator& simulator_;
    const Grid& grid_;
    bool terminal_output_;
    ModelParameters param_;
    WellModel& well_model_;
    SimulatorReportSingle failureReport_;
    std::vector<StepReport> convergence_reports_;
    ComponentName compNames_{};
    std::vector<std::vector<Scalar>> residual_norms_history_;
    Scalar current_relaxation_;
    GlobalEqVector dx_old_;
};

} // namespace Opm

#include <opm/simulators/flow/NonlinearSystem_impl.hpp>

#endif // OPM_FLOW_NONLINEAR_SYSTEM_HEADER_INCLUDED
