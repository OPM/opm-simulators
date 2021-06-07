/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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


#ifndef OPM_BLACKOILWELLMODEL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_HEADER_INCLUDED

#include <ebos/eclproblem.hh>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <cassert>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <stddef.h>

#include <opm/parser/eclipse/EclipseState/Runspec.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>

#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/WGState.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/MultisegmentWell.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/timestepping/gatherConvergenceReport.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableTerminalOutput {
    using type = UndefinedProperty;
};

} // namespace Opm::Properties

namespace Opm {

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilWellModel : public BaseAuxiliaryModule<TypeTag>
                                , public BlackoilWellModelGeneric
        {
        public:
            // ---------      Types      ---------
            typedef BlackoilModelParametersEbos<TypeTag> ModelParameters;

            using Grid = GetPropType<TypeTag, Properties::Grid>;
            using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
            using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
            using Indices = GetPropType<TypeTag, Properties::Indices>;
            using Simulator = GetPropType<TypeTag, Properties::Simulator>;
            using Scalar = GetPropType<TypeTag, Properties::Scalar>;
            using RateVector = GetPropType<TypeTag, Properties::RateVector>;
            using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
            using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

            typedef typename BaseAuxiliaryModule<TypeTag>::NeighborSet NeighborSet;

            static const int numEq = Indices::numEq;
            static const int solventSaturationIdx = Indices::solventSaturationIdx;
            static constexpr bool has_solvent_ = getPropValue<TypeTag, Properties::EnableSolvent>();
            static constexpr bool has_polymer_ = getPropValue<TypeTag, Properties::EnablePolymer>();
            static constexpr bool has_energy_ = getPropValue<TypeTag, Properties::EnableEnergy>();

            // TODO: where we should put these types, WellInterface or Well Model?
            // or there is some other strategy, like TypeTag
            typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
            typedef Dune::BlockVector<VectorBlockType> BVector;

            typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;

            typedef BlackOilPolymerModule<TypeTag> PolymerModule;

            // For the conversion between the surface volume rate and resrevoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;

            BlackoilWellModel(Simulator& ebosSimulator);

            void init();
            void initWellContainer();

            /////////////
            // <eWoms auxiliary module stuff>
            /////////////
            unsigned numDofs() const
            // No extra dofs are inserted for wells. (we use a Schur complement.)
            { return 0; }

            void addNeighbors(std::vector<NeighborSet>& neighbors) const;

            void applyInitial()
            {}

            void linearize(SparseMatrixAdapter& jacobian, GlobalEqVector& res);

            void postSolve(GlobalEqVector& deltaX)
            {
                recoverWellSolutionAndUpdateWellState(deltaX);
            }

            /////////////
            // </ eWoms auxiliary module stuff>
            /////////////

            template <class Restarter>
            void deserialize(Restarter& /* res */)
            {
                // TODO (?)
            }

            /*!
             * \brief This method writes the complete state of the well
             *        to the harddisk.
             */
            template <class Restarter>
            void serialize(Restarter& /* res*/)
            {
                // TODO (?)
            }

            void beginEpisode()
            {
                beginReportStep(ebosSimulator_.episodeIndex());
            }

            void beginTimeStep();

            void beginIteration()
            {
                assemble(ebosSimulator_.model().newtonMethod().numIterations(),
                         ebosSimulator_.timeStepSize());
            }

            void endIteration()
            { }

            void endTimeStep()
            {
                timeStepSucceeded(ebosSimulator_.time(), ebosSimulator_.timeStepSize());
            }

            void endEpisode()
            {
                endReportStep();
            }

            template <class Context>
            void computeTotalRatesForDof(RateVector& rate,
                                         const Context& context,
                                         unsigned spaceIdx,
                                         unsigned timeIdx) const;


            using WellInterfacePtr = std::shared_ptr<WellInterface<TypeTag> >;
            WellInterfacePtr well(const std::string& wellName) const;

            using BlackoilWellModelGeneric::initFromRestartFile;
            void initFromRestartFile(const RestartValue& restartValues)
            {
                initFromRestartFile(restartValues,
                                    UgGridHelpers::numCells(grid()),
                                    param_.use_multisegment_well_);
            }

            data::Wells wellData() const
            {
                auto wsrpt = this->wellState().report(UgGridHelpers::globalCell(grid()),
                                                      [this](const int well_ndex) -> bool
                                                      {
                                                          return this->wasDynamicallyShutThisTimeStep(well_ndex);
                                                      });

                this->assignWellGuideRates(wsrpt);
                this->assignShutConnections(wsrpt, this->reportStepIndex());

                return wsrpt;
            }

            // substract Binv(D)rw from r;
            void apply( BVector& r) const;

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax) const;

#if HAVE_CUDA || HAVE_OPENCL
            // accumulate the contributions of all Wells in the WellContributions object
            void getWellContributions(WellContributions& x) const;
#endif

            // apply well model with scaling of alpha
            void applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const;

            // Check if well equations is converged.
            ConvergenceReport getWellConvergence(const std::vector<Scalar>& B_avg, const bool checkGroupConvergence = false) const;

            const SimulatorReportSingle& lastReport() const;

            void addWellContributions(SparseMatrixAdapter& jacobian) const
            {
                for ( const auto& well: well_container_ ) {
                    well->addWellContributions(jacobian);
                }
            }

            // called at the beginning of a report step
            void beginReportStep(const int time_step);

            void updatePerforationIntensiveQuantities();
            // it should be able to go to prepareTimeStep(), however, the updateWellControls() and initPrimaryVariablesEvaluation()
            // makes it a little more difficult. unless we introduce if (iterationIdx != 0) to avoid doing the above functions
            // twice at the beginning of the time step
            /// Calculating the explict quantities used in the well calculation. By explicit, we mean they are cacluated
            /// at the beginning of the time step and no derivatives are included in these quantities
            void calculateExplicitQuantities(DeferredLogger& deferred_logger) const;
            // some preparation work, mostly related to group control and RESV,
            // at the beginning of each time step (Not report step)
            void prepareTimeStep(DeferredLogger& deferred_logger);
            void initPrimaryVariablesEvaluation() const;
            void updateWellControls(DeferredLogger& deferred_logger, const bool checkGroupControls);
            WellInterfacePtr getWell(const std::string& well_name) const;

        protected:
            Simulator& ebosSimulator_;

            // a vector of all the wells.
            std::vector<WellInterfacePtr > well_container_{};

            std::vector<bool> is_cell_perforated_{};

            void initializeWellState(const int           timeStepIdx,
                                     const SummaryState& summaryState);

            // create the well container
            void createWellContainer(const int time_step);

            WellInterfacePtr
            createWellPointer(const int wellID,
                              const int time_step) const;

            template <typename WellType>
            std::unique_ptr<WellType>
            createTypedWellPointer(const int wellID,
                                   const int time_step) const;

            WellInterfacePtr createWellForWellTest(const std::string& well_name, const int report_step, DeferredLogger& deferred_logger) const;


            const ModelParameters param_;
            size_t global_num_cells_{};
            // the number of the cells in the local grid
            size_t local_num_cells_{};
            double gravity_{};
            std::vector<double> depth_{};
            bool alternative_well_rate_init_{};

            std::optional<int> last_run_wellpi_{};

            std::unique_ptr<RateConverterType> rateConverter_{};

            SimulatorReportSingle last_report_{};

            // used to better efficiency of calcuation
            mutable BVector scaleAddRes_{};

            std::vector<Scalar> B_avg_{};

            const Grid& grid() const
            { return ebosSimulator_.vanguard().grid(); }

            const EclipseState& eclState() const
            { return ebosSimulator_.vanguard().eclState(); }

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble(const int iterationIdx,
                          const double dt);

            // called at the end of a time step
            void timeStepSucceeded(const double& simulationTime, const double dt);

            // called at the end of a report step
            void endReportStep();

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellState(const BVector& x);

            // setting the well_solutions_ based on well_state.
            void updatePrimaryVariables(DeferredLogger& deferred_logger);

            void setupCartesianToCompressed_(const int* global_cell, int local_num__cells);

            void updateAverageFormationFactor();

            void computePotentials(const std::size_t widx,
                                   const WellState& well_state_copy,
                                   std::string& exc_msg,
                                   ExceptionType::ExcEnum& exc_type,
                                   DeferredLogger& deferred_logger) override;

            const std::vector<double>& wellPerfEfficiencyFactors() const;

            void calculateProductivityIndexValuesShutWells(const int reportStepIdx, DeferredLogger& deferred_logger);
            void calculateProductivityIndexValues(DeferredLogger& deferred_logger);
            void calculateProductivityIndexValues(const WellInterface<TypeTag>* wellPtr,
                                                  DeferredLogger& deferred_logger);

            // The number of components in the model.
            int numComponents() const;

            int reportStepIndex() const;

            void assembleWellEq(const double dt, DeferredLogger& deferred_logger);

            void maybeDoGasLiftOptimize(DeferredLogger& deferred_logger);

            void extractLegacyCellPvtRegionIndex_();

            void extractLegacyDepth_();

            /// upate the wellTestState related to economic limits
            void updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const;

            void wellTesting(const int timeStepIdx, const double simulationTime, DeferredLogger& deferred_logger);

            void calcRates(const int fipnum,
                           const int pvtreg,
                           std::vector<double>& resv_coeff) override;

            void runWellPIScaling(const int timeStepIdx, DeferredLogger& local_deferredLogger);

             void computeWellTemperature();

        private:
            BlackoilWellModel(Simulator& ebosSimulator, const PhaseUsage& pu);
        };


} // namespace Opm

#include "BlackoilWellModel_impl.hpp"
#endif
