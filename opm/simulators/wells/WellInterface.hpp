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

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/WellState.hpp>
// NOTE: GasLiftSingleWell.hpp includes StandardWell.hpp which includes ourself
//   (WellInterface.hpp), so we need to forward declare GasLiftSingleWell
//   for it to be defined in this file. Similar for BlackoilWellModel
namespace Opm {
    template<typename TypeTag> class GasLiftSingleWell;
    template<typename TypeTag> class BlackoilWellModel;
}
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <cassert>
#include <vector>

namespace Opm
{

class WellInjectionProperties;
class WellProductionProperties;

template<typename TypeTag>
class WellInterface : public WellInterfaceIndices<GetPropType<TypeTag, Properties::FluidSystem>,
                                                  GetPropType<TypeTag, Properties::Indices>,
                                                  GetPropType<TypeTag, Properties::Scalar>>
{
    using Base = WellInterfaceIndices<GetPropType<TypeTag, Properties::FluidSystem>,
                                      GetPropType<TypeTag, Properties::Indices>,
                                      GetPropType<TypeTag, Properties::Scalar>>;
public:
    using ModelParameters = BlackoilModelParametersEbos<TypeTag>;

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using GasLiftSingleWell = ::Opm::GasLiftSingleWell<TypeTag>;
    using GLiftOptWells = typename BlackoilWellModel<TypeTag>::GLiftOptWells;
    using GLiftProdWells = typename BlackoilWellModel<TypeTag>::GLiftProdWells;
    using GLiftWellStateMap =
        typename BlackoilWellModel<TypeTag>::GLiftWellStateMap;
    using GLiftSyncGroups = typename GasLiftSingleWellGeneric::GLiftSyncGroups;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using VectorBlockType = Dune::FieldVector<Scalar, Indices::numEq>;
    using MatrixBlockType = Dune::FieldMatrix<Scalar, Indices::numEq, Indices::numEq>;
    using Eval = typename Base::Eval;
    using BVector = Dune::BlockVector<VectorBlockType>;
    using PressureMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, 1, 1>>;

    using RateConverterType =
    typename WellInterfaceFluidSystem<FluidSystem>::RateConverterType;

    using WellInterfaceFluidSystem<FluidSystem>::Gas;
    using WellInterfaceFluidSystem<FluidSystem>::Oil;
    using WellInterfaceFluidSystem<FluidSystem>::Water;

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
                  const ParallelWellInfo& pw_info,
                  const int time_step,
                  const ModelParameters& param,
                  const RateConverterType& rate_converter,
                  const int pvtRegionIdx,
                  const int num_components,
                  const int num_phases,
                  const int index_of_well,
                  const std::vector<PerforationData>& perf_data);

    /// Virtual destructor
    virtual ~WellInterface() = default;

    virtual void init(const PhaseUsage* phase_usage_arg,
                      const std::vector<double>& depth_arg,
                      const double gravity_arg,
                      const int num_cells,
                      const std::vector< Scalar >& B_avg,
                      const bool changed_to_open_this_step);

    virtual void initPrimaryVariablesEvaluation() = 0;

    virtual ConvergenceReport getWellConvergence(const SummaryState& summary_state,
                                                 const WellState& well_state,
                                                 const std::vector<double>& B_avg,
                                                 DeferredLogger& deferred_logger,
                                                 const bool relax_tolerance) const = 0;

    virtual void solveEqAndUpdateWellState(const SummaryState& summary_state,
                                           WellState& well_state,
                                           DeferredLogger& deferred_logger) = 0;

    void assembleWellEq(const Simulator& ebosSimulator,
                        const double dt,
                        WellState& well_state,
                        const GroupState& group_state,
                        DeferredLogger& deferred_logger);

    void assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                        const double dt,
                                        WellState& well_state,
                                        const GroupState& group_state,
                                        DeferredLogger& deferred_logger);

    // TODO: better name or further refactoring the function to make it more clear
    void prepareWellBeforeAssembling(const Simulator& ebosSimulator,
                                     const double dt,
                                     WellState& well_state,
                                     const GroupState& group_state,
                                     DeferredLogger& deferred_logger);


    virtual void computeWellRatesWithBhp(
        const Simulator& ebosSimulator,
        const double& bhp,
        std::vector<double>& well_flux,
        DeferredLogger& deferred_logger
    ) const = 0;

    virtual std::optional<double> computeBhpAtThpLimitProdWithAlq(
        const Simulator& ebos_simulator,
        const SummaryState& summary_state,
        const double alq_value,
        DeferredLogger& deferred_logger
    ) const = 0;

    /// using the solution x to recover the solution xw for wells and applying
    /// xw to update Well State
    virtual void recoverWellSolutionAndUpdateWellState(const SummaryState& summary_state,
                                                       const BVector& x,
                                                       WellState& well_state,
                                                       DeferredLogger& deferred_logger) = 0;

    /// Ax = Ax - C D^-1 B x
    virtual void apply(const BVector& x, BVector& Ax) const = 0;

    /// r = r - C D^-1 Rw
    virtual void apply(BVector& r) const = 0;

    // TODO: before we decide to put more information under mutable, this function is not const
    virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                       const WellState& well_state,
                                       std::vector<double>& well_potentials,
                                       DeferredLogger& deferred_logger) = 0;

    virtual void updateWellStateWithTarget(const Simulator& ebos_simulator,
                                           const GroupState& group_state,
                                           WellState& well_state,
                                           DeferredLogger& deferred_logger) const;

    virtual void computeWellRatesWithBhpIterations(const Simulator& ebosSimulator,
                                                   const Scalar& bhp,
                                                   std::vector<double>& well_flux,
                                                   DeferredLogger& deferred_logger) const = 0;

    bool updateWellStateWithTHPTargetProd(const Simulator& ebos_simulator,
                                          WellState& well_state,
                                          DeferredLogger& deferred_logger) const;

    enum class IndividualOrGroup { Individual, Group, Both };
    bool updateWellControl(const Simulator& ebos_simulator,
                           const IndividualOrGroup iog,
                           WellState& well_state,
                           const GroupState& group_state,
                           DeferredLogger& deferred_logger) /* const */;

    bool updateWellControlAndStatusLocalIteration(const Simulator& ebos_simulator,
                                                  WellState& well_state,
                                                  const GroupState& group_state,
                                                  const Well::InjectionControls& inj_controls,
                                                  const Well::ProductionControls& prod_controls,
                                                  const double WQTotal,
                                                  DeferredLogger& deferred_logger);

    virtual void updatePrimaryVariables(const SummaryState& summary_state,
                                        const WellState& well_state,
                                        DeferredLogger& deferred_logger) = 0;

    virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                             const WellState& well_state,
                                             DeferredLogger& deferred_logger) = 0; // should be const?

    virtual void updateProductivityIndex(const Simulator& ebosSimulator,
                                         const WellProdIndexCalculator& wellPICalc,
                                         WellState& well_state,
                                         DeferredLogger& deferred_logger) const = 0;

    virtual double connectionDensity(const int globalConnIdx,
                                     const int openConnIdx) const = 0;

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
                                          const WellState& well_state) const = 0;

    void addCellRates(RateVector& rates, int cellIdx) const;

    Scalar volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const;

    // TODO: theoretically, it should be a const function
    // Simulator is not const is because that assembleWellEq is non-const Simulator
    void wellTesting(const Simulator& simulator,
                     const double simulation_time,
                     /* const */ WellState& well_state, const GroupState& group_state, WellTestState& welltest_state,
                     DeferredLogger& deferred_logger);

    void checkWellOperability(const Simulator& ebos_simulator, const WellState& well_state, DeferredLogger& deferred_logger);

    bool gliftBeginTimeStepWellTestIterateWellEquations(
        const Simulator& ebos_simulator,
        const double dt,
        WellState& well_state,
        const GroupState &group_state,
        DeferredLogger& deferred_logger);

    void gliftBeginTimeStepWellTestUpdateALQ(const Simulator& ebos_simulator,
                                             WellState& well_state,
                                             DeferredLogger& deferred_logger);

    // check whether the well is operable under the current reservoir condition
    // mostly related to BHP limit and THP limit
    void updateWellOperability(const Simulator& ebos_simulator,
                               const WellState& well_state,
                               DeferredLogger& deferred_logger);

    bool updateWellOperabilityFromWellEq(const Simulator& ebos_simulator,
                                         const WellState& well_state,
                                         DeferredLogger& deferred_logger);

    // update perforation water throughput based on solved water rate
    virtual void updateWaterThroughput(const double dt, WellState& well_state) const = 0;

    /// Compute well rates based on current reservoir conditions and well variables.
    /// Used in updateWellStateRates().
    virtual std::vector<double> computeCurrentWellRates(const Simulator& ebosSimulator,
                                                        DeferredLogger& deferred_logger) const = 0;

    /// Modify the well_state's rates if there is only one nonzero rate.
    /// If so, that rate is kept as is, but the others are set proportionally
    /// to the rates returned by computeCurrentWellRates().
    void updateWellStateRates(const Simulator& ebosSimulator,
                              WellState& well_state,
                              DeferredLogger& deferred_logger) const;

    void solveWellEquation(const Simulator& ebosSimulator,
                           WellState& well_state,
                           const GroupState& group_state,
                           DeferredLogger& deferred_logger);

    const std::vector<RateVector>& connectionRates() const
    {
        return connectionRates_;
    }

    virtual std::vector<double> getPrimaryVars() const
    {
        return {};
    }

    virtual int setPrimaryVars(std::vector<double>::const_iterator)
    {
        return 0;
    }

    std::vector<double> wellIndex(const int perf, const IntensiveQuantities& intQuants, const double trans_mult, const SingleWellState& ws) const;

    void updateConnectionDFactor(const Simulator& simulator, SingleWellState& ws) const;

    void updateConnectionTransmissibilityFactor(const Simulator& simulator, SingleWellState& ws) const;


protected:
    // simulation parameters
    const ModelParameters& param_;
    std::vector<RateVector> connectionRates_;
    std::vector< Scalar > B_avg_;
    bool changed_to_stopped_this_step_ = false;
    bool thp_update_iterations = false;

    double wpolymer() const;

    double wfoam() const;

    double wsalt() const;

    double wmicrobes() const;

    double woxygen() const;

    double wurea() const;

    virtual double getRefDensity() const = 0;

    // Component fractions for each phase for the well
    const std::vector<double>& compFrac() const;

    std::vector<double> initialWellRateFractions(const Simulator& ebosSimulator, const WellState& well_state) const;

    // check whether the well is operable under BHP limit with current reservoir condition
    virtual void checkOperabilityUnderBHPLimit(const WellState& well_state, const Simulator& ebos_simulator, DeferredLogger& deferred_logger) =0;

    // check whether the well is operable under THP limit with current reservoir condition
    virtual void checkOperabilityUnderTHPLimit(const Simulator& ebos_simulator, const WellState& well_state, DeferredLogger& deferred_logger) =0;

    virtual void updateIPR(const Simulator& ebos_simulator, DeferredLogger& deferred_logger) const=0;

    virtual void assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                                const double dt,
                                                const WellInjectionControls& inj_controls,
                                                const WellProductionControls& prod_controls,
                                                WellState& well_state,
                                                const GroupState& group_state,
                                                DeferredLogger& deferred_logger) = 0;

    // iterate well equations with the specified control until converged
    virtual bool iterateWellEqWithControl(const Simulator& ebosSimulator,
                                          const double dt,
                                          const WellInjectionControls& inj_controls,
                                          const WellProductionControls& prod_controls,
                                          WellState& well_state,
                                          const GroupState& group_state,
                                          DeferredLogger& deferred_logger) = 0;

    virtual bool iterateWellEqWithSwitching(const Simulator& ebosSimulator,
                                            const double dt,
                                            const WellInjectionControls& inj_controls,
                                            const WellProductionControls& prod_controls,
                                            WellState& well_state,
                                            const GroupState& group_state,
                                            DeferredLogger& deferred_logger) = 0;

    bool iterateWellEquations(const Simulator& ebosSimulator,
                              const double dt,
                              WellState& well_state,
                              const GroupState& group_state,
                              DeferredLogger& deferred_logger);

    bool solveWellForTesting(const Simulator& ebosSimulator, WellState& well_state, const GroupState& group_state,
                             DeferredLogger& deferred_logger);

    Eval getPerfCellPressure(const FluidState& fs) const;

    // get the mobility for specific perforation
    template<class Value, class Callback>
    void getMobility(const Simulator& ebosSimulator,
                     const int perf,
                     std::vector<Value>& mob,
                     Callback& extendEval,
                     [[maybe_unused]] DeferredLogger& deferred_logger) const;

    void computeConnLevelProdInd(const FluidState& fs,
                                 const std::function<double(const double)>& connPICalc,
                                 const std::vector<Scalar>& mobility,
                                 double* connPI) const;

    void computeConnLevelInjInd(const FluidState& fs,
                                const Phase preferred_phase,
                                const std::function<double(const double)>& connIICalc,
                                const std::vector<Scalar>& mobility,
                                double* connII,
                                DeferredLogger& deferred_logger) const;

    double computeConnectionDFactor(const int perf, const IntensiveQuantities& intQuants, const double trans_mult, const double total_tw, const SingleWellState& ws) const;


};

}

#include "WellInterface_impl.hpp"

#endif // OPM_WELLINTERFACE_HEADER_INCLUDED
