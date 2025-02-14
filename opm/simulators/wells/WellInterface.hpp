/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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

#ifndef OPM_WELLINTERFACE_HEADER_INCLUDED
#define OPM_WELLINTERFACE_HEADER_INCLUDED

// NOTE: GasLiftSingleWell.hpp includes StandardWell.hpp which includes ourself
//   (WellInterface.hpp), so we need to forward declare GasLiftSingleWell
//   for it to be defined in this file. Similar for BlackoilWellModel
namespace Opm {
    template<typename TypeTag> class GasLiftSingleWell;
    template<typename TypeTag> class BlackoilWellModel;
}

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/material/fluidstates/BlackOilFluidState.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/simulators/linalg/linalgproperties.hh>

#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <opm/simulators/utils/BlackoilPhases.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Evaluation.hpp>

#include <vector>

namespace Opm
{

class WellInjectionProperties;
class WellProductionProperties;

template<typename TypeTag>
class WellInterface : public WellInterfaceIndices<GetPropType<TypeTag, Properties::FluidSystem>,
                                                  GetPropType<TypeTag, Properties::Indices>>
{
    using Base = WellInterfaceIndices<GetPropType<TypeTag, Properties::FluidSystem>,
                                      GetPropType<TypeTag, Properties::Indices>>;
public:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using GasLiftSingleWell = ::Opm::GasLiftSingleWell<TypeTag>;
    using GLiftEclWells = typename GasLiftGroupInfo<Scalar>::GLiftEclWells;

    using VectorBlockType = Dune::FieldVector<Scalar, Indices::numEq>;
    using MatrixBlockType = Dune::FieldMatrix<Scalar, Indices::numEq, Indices::numEq>;
    using Eval = typename Base::Eval;
    using BVector = Dune::BlockVector<VectorBlockType>;
    using PressureMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<Scalar, 1, 1>>;

    using RateConverterType =
    typename WellInterfaceFluidSystem<FluidSystem>::RateConverterType;

    using WellInterfaceFluidSystem<FluidSystem>::Gas;
    using WellInterfaceFluidSystem<FluidSystem>::Oil;
    using WellInterfaceFluidSystem<FluidSystem>::Water;

    using ModelParameters = typename Base::ModelParameters;

    static constexpr bool has_solvent = getPropValue<TypeTag, Properties::EnableSolvent>();
    static constexpr bool has_zFraction = getPropValue<TypeTag, Properties::EnableExtbo>();
    static constexpr bool has_polymer = getPropValue<TypeTag, Properties::EnablePolymer>();
    static constexpr bool has_energy = getPropValue<TypeTag, Properties::EnableEnergy>();
    static const bool has_temperature = getPropValue<TypeTag, Properties::EnableTemperature>();
    // flag for polymer molecular weight related
    static constexpr bool has_polymermw = getPropValue<TypeTag, Properties::EnablePolymerMW>();
    static constexpr bool has_foam = getPropValue<TypeTag, Properties::EnableFoam>();
    static constexpr bool has_brine = getPropValue<TypeTag, Properties::EnableBrine>();
    static constexpr bool has_watVapor = getPropValue<TypeTag, Properties::EnableVapwat>();
    static constexpr bool has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>();
    static constexpr bool has_saltPrecip = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();
    static constexpr bool has_micp = getPropValue<TypeTag, Properties::EnableMICP>();

    // For the conversion between the surface volume rate and reservoir voidage rate
    using FluidState = BlackOilFluidState<Eval,
                                          FluidSystem,
                                          has_temperature,
                                          has_energy,
                                          Indices::compositionSwitchIdx >= 0,
                                          has_watVapor,
                                          has_brine,
                                          has_saltPrecip,
                                          has_disgas_in_water,
                                          Indices::numPhases >;
    /// Constructor
    WellInterface(const Well& well,
                  const ParallelWellInfo<Scalar>& pw_info,
                  const int time_step,
                  const ModelParameters& param,
                  const RateConverterType& rate_converter,
                  const int pvtRegionIdx,
                  const int num_components,
                  const int num_phases,
                  const int index_of_well,
                  const std::vector<PerforationData<Scalar>>& perf_data);

    /// Virtual destructor
    virtual ~WellInterface() = default;

    virtual void init(const PhaseUsage* phase_usage_arg,
                      const std::vector<Scalar>& depth_arg,
                      const Scalar gravity_arg,
                      const std::vector<Scalar>& B_avg,
                      const bool changed_to_open_this_step);

    virtual ConvergenceReport getWellConvergence(const Simulator& simulator,
                                                 const WellState<Scalar>& well_state,
                                                 const std::vector<Scalar>& B_avg,
                                                 DeferredLogger& deferred_logger,
                                                 const bool relax_tolerance) const = 0;

    virtual void solveEqAndUpdateWellState(const Simulator& simulator,
                                           WellState<Scalar>& well_state,
                                           DeferredLogger& deferred_logger) = 0;

    void assembleWellEq(const Simulator& simulator,
                        const double dt,
                        WellState<Scalar>& well_state,
                        const GroupState<Scalar>& group_state,
                        DeferredLogger& deferred_logger);

    void assembleWellEqWithoutIteration(const Simulator& simulator,
                                        const double dt,
                                        WellState<Scalar>& well_state,
                                        const GroupState<Scalar>& group_state,
                                        DeferredLogger& deferred_logger);

    // TODO: better name or further refactoring the function to make it more clear
    void prepareWellBeforeAssembling(const Simulator& simulator,
                                     const double dt,
                                     WellState<Scalar>& well_state,
                                     const GroupState<Scalar>& group_state,
                                     DeferredLogger& deferred_logger);


    virtual void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                         const Scalar& bhp,
                                         std::vector<Scalar>& well_flux,
                                         DeferredLogger& deferred_logger) const = 0;

    virtual std::optional<Scalar>
    computeBhpAtThpLimitProdWithAlq(const Simulator& ebos_simulator,
                                    const SummaryState& summary_state,
                                    const Scalar alq_value,
                                    DeferredLogger& deferred_logger,
                                    bool iterate_if_no_solution) const = 0;

    /// using the solution x to recover the solution xw for wells and applying
    /// xw to update Well State
    virtual void recoverWellSolutionAndUpdateWellState(const Simulator& simulator,
                                                       const BVector& x,
                                                       WellState<Scalar>& well_state,
                                                       DeferredLogger& deferred_logger) = 0;

    /// Ax = Ax - C D^-1 B x
    virtual void apply(const BVector& x, BVector& Ax) const = 0;

    /// r = r - C D^-1 Rw
    virtual void apply(BVector& r) const = 0;

    // TODO: before we decide to put more information under mutable, this function is not const
    virtual void computeWellPotentials(const Simulator& simulator,
                                       const WellState<Scalar>& well_state,
                                       std::vector<Scalar>& well_potentials,
                                       DeferredLogger& deferred_logger) = 0;

    virtual void updateWellStateWithTarget(const Simulator& simulator,
                                           const GroupState<Scalar>& group_state,
                                           WellState<Scalar>& well_state,
                                           DeferredLogger& deferred_logger) const;

    virtual void computeWellRatesWithBhpIterations(const Simulator& simulator,
                                                   const Scalar& bhp,
                                                   std::vector<Scalar>& well_flux,
                                                   DeferredLogger& deferred_logger) const = 0;

    bool wellUnderZeroRateTarget(const Simulator& simulator,
                                 const WellState<Scalar>& well_state,
                                 DeferredLogger& deferred_logger) const;

    bool wellUnderZeroGroupRateTarget(const Simulator& simulator,
                                      const WellState<Scalar>& well_state,
                                      DeferredLogger& deferred_logger,
                                      std::optional<bool> group_control = std::nullopt) const;

    bool stoppedOrZeroRateTarget(const Simulator& simulator,
                                 const WellState<Scalar>& well_state,
                                 DeferredLogger& deferred_logger) const;

    bool updateWellStateWithTHPTargetProd(const Simulator& simulator,
                                          WellState<Scalar>& well_state,
                                          DeferredLogger& deferred_logger) const;

    enum class IndividualOrGroup { Individual, Group, Both };
    bool updateWellControl(const Simulator& simulator,
                           const IndividualOrGroup iog,
                           WellState<Scalar>& well_state,
                           const GroupState<Scalar>& group_state,
                           DeferredLogger& deferred_logger) /* const */;

    bool updateWellControlAndStatusLocalIteration(const Simulator& simulator,
                                                  WellState<Scalar>& well_state,
                                                  const GroupState<Scalar>& group_state,
                                                  const Well::InjectionControls& inj_controls,
                                                  const Well::ProductionControls& prod_controls,
                                                  const Scalar WQTotal,
                                                  DeferredLogger& deferred_logger, 
                                                  const bool fixed_control = false, 
                                                  const bool fixed_status = false);

    virtual void updatePrimaryVariables(const Simulator& simulator,
                                        const WellState<Scalar>& well_state,
                                        DeferredLogger& deferred_logger) = 0;

    virtual void calculateExplicitQuantities(const Simulator& simulator,
                                             const WellState<Scalar>& well_state,
                                             DeferredLogger& deferred_logger) = 0; // should be const?

    virtual void updateProductivityIndex(const Simulator& simulator,
                                         const WellProdIndexCalculator<Scalar>& wellPICalc,
                                         WellState<Scalar>& well_state,
                                         DeferredLogger& deferred_logger) const = 0;

    /// \brief Wether the Jacobian will also have well contributions in it.
    virtual bool jacobianContainsWellContributions() const
    {
        return false;
    }

    // Add well contributions to matrix
    virtual void addWellContributions(SparseMatrixAdapter&) const = 0;

    virtual void addWellPressureEquations(PressureMatrix& mat,
                                          const BVector& x,
                                          const int pressureVarIndex,
                                          const bool use_well_weights,
                                          const WellState<Scalar>& well_state) const = 0;

    void addCellRates(RateVector& rates, int cellIdx) const;

    Scalar volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const;

    // TODO: theoretically, it should be a const function
    // Simulator is not const is because that assembleWellEq is non-const Simulator
    void wellTesting(const Simulator& simulator,
                     const double simulation_time,
                     /* const */ WellState<Scalar>& well_state,
                     const GroupState<Scalar>& group_state,
                     WellTestState& welltest_state,
                     const PhaseUsage& phase_usage,
                     GLiftEclWells& ecl_well_map,
                     std::map<std::string, double>& open_times,
                     DeferredLogger& deferred_logger);

    void checkWellOperability(const Simulator& simulator,
                              const WellState<Scalar>& well_state,
                              DeferredLogger& deferred_logger);

    bool gliftBeginTimeStepWellTestIterateWellEquations(const Simulator& ebos_simulator,
                                                        const double dt,
                                                        WellState<Scalar>& well_state,
                                                        const GroupState<Scalar>& group_state,
                                                        DeferredLogger& deferred_logger);

    void gliftBeginTimeStepWellTestUpdateALQ(const Simulator& simulator,
                                             WellState<Scalar>& well_state,
                                             const GroupState<Scalar>& group_state,
                                             const PhaseUsage& phase_usage,
                                             GLiftEclWells& ecl_well_map,
                                             DeferredLogger& deferred_logger);

    // check whether the well is operable under the current reservoir condition
    // mostly related to BHP limit and THP limit
    void updateWellOperability(const Simulator& simulator,
                               const WellState<Scalar>& well_state,
                               DeferredLogger& deferred_logger);

    bool updateWellOperabilityFromWellEq(const Simulator& simulator,
                                         const WellState<Scalar>& well_state,
                                         DeferredLogger& deferred_logger);

    // update perforation water throughput based on solved water rate
    virtual void updateWaterThroughput(const double dt,
                                       WellState<Scalar>& well_state) const = 0;

    /// Compute well rates based on current reservoir conditions and well variables.
    /// Used in updateWellStateRates().
    virtual std::vector<Scalar>
    computeCurrentWellRates(const Simulator& simulator,
                            DeferredLogger& deferred_logger) const = 0;

    /// Modify the well_state's rates if there is only one nonzero rate.
    /// If so, that rate is kept as is, but the others are set proportionally
    /// to the rates returned by computeCurrentWellRates().
    void initializeWellState(const Simulator& simulator,
                             const GroupState<Scalar>& group_state,
                             const SummaryState& summary_state,
                             WellState<Scalar>& well_state,
                             DeferredLogger& deferred_logger) const;

    // void updateWellStateRates(const Simulator& simulator,
    //                           WellState<Scalar>& well_state,
    //                           DeferredLogger& deferred_logger) const;

    void solveWellEquation(const Simulator& simulator,
                           WellState<Scalar>& well_state,
                           const GroupState<Scalar>& group_state,
                           DeferredLogger& deferred_logger);

    const std::vector<RateVector>& connectionRates() const
    {
        return connectionRates_;
    }

    std::vector<Scalar> wellIndex(const int perf,
                                  const IntensiveQuantities& intQuants,
                                  const Scalar trans_mult,
                                  const SingleWellState<Scalar>& ws) const;

    void updateConnectionDFactor(const Simulator& simulator,
                                 SingleWellState<Scalar>& ws) const;

    void updateConnectionTransmissibilityFactor(const Simulator& simulator,
                                                SingleWellState<Scalar>& ws) const;

    virtual bool iterateWellEqWithSwitching(const Simulator& simulator,
                                            const double dt,
                                            const WellInjectionControls& inj_controls,
                                            const WellProductionControls& prod_controls,
                                            WellState<Scalar>& well_state,
                                            const GroupState<Scalar>& group_state,
                                            DeferredLogger& deferred_logger, 
                                            const bool fixed_control = false, 
                                            const bool fixed_status = false) = 0;
protected:
    // simulation parameters
    std::vector<RateVector> connectionRates_;
    std::vector<Scalar> B_avg_;
    bool changed_to_stopped_this_step_ = false;
    bool thp_update_iterations = false;

    Scalar wpolymer() const;
    Scalar wfoam() const;
    Scalar wsalt() const;
    Scalar wmicrobes() const;
    Scalar woxygen() const;
    Scalar wurea() const;

    virtual Scalar getRefDensity() const = 0;

    // Component fractions for each phase for the well
    const std::vector<Scalar>& compFrac() const;

    std::vector<Scalar>
    initialWellRateFractions(const Simulator& ebosSimulator,
                             const WellState<Scalar>& well_state) const;

    // check whether the well is operable under BHP limit with current reservoir condition
    virtual void checkOperabilityUnderBHPLimit(const WellState<Scalar>& well_state,
                                               const Simulator& simulator,
                                               DeferredLogger& deferred_logger) = 0;

    // check whether the well is operable under THP limit with current reservoir condition
    virtual void checkOperabilityUnderTHPLimit(const Simulator& simulator,
                                               const WellState<Scalar>& well_state,
                                               DeferredLogger& deferred_logger) = 0;

    virtual void updateIPR(const Simulator& simulator,
                           DeferredLogger& deferred_logger) const = 0;

    virtual void assembleWellEqWithoutIteration(const Simulator& simulator,
                                                const double dt,
                                                const WellInjectionControls& inj_controls,
                                                const WellProductionControls& prod_controls,
                                                WellState<Scalar>& well_state,
                                                const GroupState<Scalar>& group_state,
                                                DeferredLogger& deferred_logger) = 0;

    // iterate well equations with the specified control until converged
    virtual bool iterateWellEqWithControl(const Simulator& simulator,
                                          const double dt,
                                          const WellInjectionControls& inj_controls,
                                          const WellProductionControls& prod_controls,
                                          WellState<Scalar>& well_state,
                                          const GroupState<Scalar>& group_state,
                                          DeferredLogger& deferred_logger) = 0;

    virtual void updateIPRImplicit(const Simulator& simulator,
                                   WellState<Scalar>& well_state,
                                   DeferredLogger& deferred_logger) = 0;                                            

    bool iterateWellEquations(const Simulator& simulator,
                              const double dt,
                              WellState<Scalar>& well_state,
                              const GroupState<Scalar>& group_state,
                              DeferredLogger& deferred_logger);

    bool solveWellWithTHPConstraint(const Simulator& simulator,
                                    const double dt,
                                    const Well::InjectionControls& inj_controls,
                                    const Well::ProductionControls& prod_controls,
                                    WellState<Scalar>& well_state,
                                    const GroupState<Scalar>& group_state,
                                    DeferredLogger& deferred_logger);

    std::optional<Scalar>
    estimateOperableBhp(const Simulator& ebos_simulator,
                        const double dt,
                        WellState<Scalar>& well_state,
                        const SummaryState& summary_state,
                        DeferredLogger& deferred_logger);

    bool solveWellWithBhp(const Simulator& simulator,
                          const double dt,
                          const Scalar bhp,
                          WellState<Scalar>& well_state,
                          DeferredLogger& deferred_logger);         

    bool solveWellWithZeroRate(const Simulator& simulator,
                               const double dt,
                               WellState<Scalar>& well_state,
                               DeferredLogger& deferred_logger);                                                                                                       

    bool solveWellForTesting(const Simulator& simulator,
                             WellState<Scalar>& well_state,
                             const GroupState<Scalar>& group_state,
                             DeferredLogger& deferred_logger);
    

    template<class GasLiftSingleWell>
    std::unique_ptr<GasLiftSingleWell> initializeGliftWellTest_(const Simulator& simulator,
                                                                WellState<Scalar>& well_state,
                                                                const GroupState<Scalar>& group_state,
                                                                const PhaseUsage& phase_usage,
                                                                GLiftEclWells& ecl_well_map,
                                                                DeferredLogger& deferred_logger);

    Eval getPerfCellPressure(const FluidState& fs) const;

    // get the mobility for specific perforation
    template<class Value, class Callback>
    void getMobility(const Simulator& simulator,
                     const int perf,
                     std::vector<Value>& mob,
                     Callback& extendEval,
                     [[maybe_unused]] DeferredLogger& deferred_logger) const;

    void computeConnLevelProdInd(const FluidState& fs,
                                 const std::function<Scalar(const Scalar)>& connPICalc,
                                 const std::vector<Scalar>& mobility,
                                 Scalar* connPI) const;

    void computeConnLevelInjInd(const FluidState& fs,
                                const Phase preferred_phase,
                                const std::function<Scalar(const Scalar)>& connIICalc,
                                const std::vector<Scalar>& mobility,
                                Scalar* connII,
                                DeferredLogger& deferred_logger) const;

    Scalar computeConnectionDFactor(const int perf,
                                    const IntensiveQuantities& intQuants,
                                    const SingleWellState<Scalar>& ws) const;
};

} // namespace Opm

#include "WellInterface_impl.hpp"

#endif // OPM_WELLINTERFACE_HEADER_INCLUDED
