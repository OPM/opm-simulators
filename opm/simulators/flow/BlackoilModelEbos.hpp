/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#ifndef OPM_BLACKOILMODELEBOS_HEADER_INCLUDED
#define OPM_BLACKOILMODELEBOS_HEADER_INCLUDED

#include <ebos/eclproblem.hh>
#include <opm/models/utils/start.hh>

#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>

#include <opm/simulators/flow/NonlinearSolverEbos.hpp>
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/wells/WellConnectionAuxiliaryModule.hpp>
#include <opm/simulators/flow/aspinPartition.hpp>
#include <opm/simulators/flow/SubDomain.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/linalg/extractMatrix.hpp>

#include <opm/grid/common/SubGridView.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/grid/UnstructuredGrid.h>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/flow/NonlinearSolverEbos.hpp>
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>
#include <opm/simulators/linalg/ISTLSolverEbos.hpp>



#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/WellConnectionAuxiliaryModule.hpp>


#include <dune/istl/matrixmarket.hh>
#include <opm/simulators/linalg/MatrixMarketSpecializations.hpp>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/timer.hh>
#include <dune/common/unused.hh>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace Opm::Properties {

namespace TTag {
struct EclFlowProblem {
    using InheritsFrom = std::tuple<FlowTimeSteppingParameters, FlowModelParameters,
                                    FlowNonLinearSolver, EclBaseProblem, BlackOilModel>;
};
}
template<class TypeTag>
struct OutputDir<TypeTag, TTag::EclFlowProblem> {
    static constexpr auto value = "";
};
template<class TypeTag>
struct EnableDebuggingChecks<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};
// default in flow is to formulate the equations in surface volumes
template<class TypeTag>
struct BlackoilConserveSurfaceVolume<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EclAquiferModel<TypeTag, TTag::EclFlowProblem> {
    using type = BlackoilAquiferModel<TypeTag>;
};

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableTemperature<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableFoam<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableBrine<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableSaltPrecipitation<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableMICP<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EclWellModel<TypeTag, TTag::EclFlowProblem> {
    using type = BlackoilWellModel<TypeTag>;
};
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::EclFlowProblem> {
    using type = TTag::FlowIstlSolver;
};

} // namespace Opm::Properties

namespace Opm {
    /// A model implementation for three-phase black oil.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa. It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    template <class TypeTag>
    class BlackoilModelEbos
    {
    public:
        // ---------  Types and enums  ---------
        typedef BlackoilModelParametersEbos<TypeTag> ModelParameters;

        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
        using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

        typedef double Scalar;
        static const int numEq = Indices::numEq;
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiZfracEqIdx = Indices::contiZfracEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        static const int contiEnergyEqIdx = Indices::contiEnergyEqIdx;
        static const int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;
        static const int contiFoamEqIdx = Indices::contiFoamEqIdx;
	static const int contiBrineEqIdx = Indices::contiBrineEqIdx;
        static const int contiMicrobialEqIdx = Indices::contiMicrobialEqIdx;
        static const int contiOxygenEqIdx = Indices::contiOxygenEqIdx;
        static const int contiUreaEqIdx = Indices::contiUreaEqIdx;
        static const int contiBiofilmEqIdx = Indices::contiBiofilmEqIdx;
        static const int contiCalciteEqIdx = Indices::contiCalciteEqIdx;
        static const int solventSaturationIdx = Indices::solventSaturationIdx;
        static const int zFractionIdx = Indices::zFractionIdx;
        static const int polymerConcentrationIdx = Indices::polymerConcentrationIdx;
        static const int polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;
        static const int temperatureIdx = Indices::temperatureIdx;
        static const int foamConcentrationIdx = Indices::foamConcentrationIdx;
	static const int saltConcentrationIdx = Indices::saltConcentrationIdx;
        static const int microbialConcentrationIdx = Indices::microbialConcentrationIdx;
        static const int oxygenConcentrationIdx = Indices::oxygenConcentrationIdx;
        static const int ureaConcentrationIdx = Indices::ureaConcentrationIdx;
        static const int biofilmConcentrationIdx = Indices::biofilmConcentrationIdx;
        static const int calciteConcentrationIdx = Indices::calciteConcentrationIdx;

        typedef Dune::FieldVector<Scalar, numEq >        VectorBlockType;
        typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
        typedef typename SparseMatrixAdapter::IstlMatrix Mat;
        typedef Dune::BlockVector<VectorBlockType>      BVector;

        using Domain = SubDomain<Grid>;

        typedef ISTLSolverEbos<TypeTag> ISTLSolverType;

        class ComponentName
        {
        public:
            ComponentName()
                : names_(numEq)
            {
                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }

                    const unsigned canonicalCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
                    names_[Indices::canonicalToActiveComponentIndex(canonicalCompIdx)]
                        = FluidSystem::componentName(canonicalCompIdx);
                }

                if constexpr (has_solvent_) {
                    names_[solventSaturationIdx] = "Solvent";
                }

                if constexpr (has_extbo_) {
                    names_[zFractionIdx] = "ZFraction";
                }

                if constexpr (has_polymer_) {
                    names_[polymerConcentrationIdx] = "Polymer";
                }

                if constexpr (has_polymermw_) {
                    assert(has_polymer_);
                    names_[polymerMoleWeightIdx] = "MolecularWeightP";
                }

                if constexpr (has_energy_) {
                    names_[temperatureIdx] = "Energy";
                }

                if constexpr (has_foam_) {
                    names_[foamConcentrationIdx] = "Foam";
                }

                if constexpr (has_brine_) {
                    names_[saltConcentrationIdx] = "Brine";
                }

                if constexpr (has_micp_) {
                    names_[microbialConcentrationIdx] = "Microbes";
                    names_[oxygenConcentrationIdx] = "Oxygen";
                    names_[ureaConcentrationIdx] = "Urea";
                    names_[biofilmConcentrationIdx] = "Biofilm";
                    names_[calciteConcentrationIdx] = "Calcite";
                }
            }

            const std::string& name(const int compIdx) const
            {
                return this->names_[compIdx];
            }

        private:
            std::vector<std::string> names_{};
        };

        //typedef typename SolutionVector :: value_type            PrimaryVariables ;

        // ---------  Public methods  ---------

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] wells            well structure
        /// \param[in] vfp_properties   Vertical flow performance tables
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilModelEbos(Simulator& ebosSimulator,
                          const ModelParameters& param,
                          BlackoilWellModel<TypeTag>& well_model,
                          const bool terminal_output)
        : ebosSimulator_(ebosSimulator)
        , grid_(ebosSimulator_.vanguard().grid())
        , phaseUsage_(phaseUsageFromDeck(eclState()))
        , param_( param )
        , well_model_ (well_model)
        , terminal_output_ (terminal_output)
        , current_relaxation_(1.0)
        , dx_old_(ebosSimulator_.model().numGridDof())
        {
            // compute global sum of number of cells
            global_nc_ = detail::countGlobalCells(grid_);
            convergence_reports_.reserve(300); // Often insufficient, but avoids frequent moves.
            // TODO: remember to fix!
            if (param_.nonlinear_solver_ == "aspin" || param_.nonlinear_solver_ == "nldd" || param_.nonlinear_solver_ == "purelocal") {
                setupAspinDomains();
            }
        }


        bool isParallel() const
        { return  grid_.comm().size() > 1; }


        const EclipseState& eclState() const
        { return ebosSimulator_.vanguard().eclState(); }



        void setupAspinDomains()
        {
            // Create partitions.
            const int num_cells = detail::countLocalInteriorCells(grid_);
            auto [partition_vector, num_partitions] = partitionCells(num_cells);

            // Scan through partitioning to get correct size for each.
            std::vector<int> sizes(num_partitions, 0);
            for (int p : partition_vector) {
                ++sizes[p];
            }

            // Set up correctly sized vectors of entity seeds and of indices for each partition.
            using EntitySeed = typename Grid::template Codim<0>::EntitySeed;
            std::vector<std::vector<EntitySeed>> seeds(num_partitions);
            std::vector<std::vector<int>> partitions(num_partitions);
            for (int part = 0; part < num_partitions; ++part) {
                seeds[part].resize(sizes[part]);
                partitions[part].resize(sizes[part]);
            }

            // Iterate through grid once, setting the seeds of all partitions.
            std::vector<int> count(num_partitions, 0);
            const auto beg = grid_.template leafbegin<0>();
            const auto end = grid_.template leafend<0>();
            int cell = 0;
            for (auto it = beg; it != end; ++it, ++cell) {
                const int p = partition_vector[cell];
                seeds[p][count[p]] = it->seed();
                partitions[p][count[p]] = cell;
                ++count[p];
            }
            assert(cell == num_cells);
            assert(count == sizes);

            // Create the domains.
            const int num_domains = partitions.size();
            for (int index = 0; index < num_domains; ++index) {
                std::vector<bool> interior(num_cells, false);
                for (int ix : partitions[index]) {
                    interior[ix] = true;
                }
                Dune::SubGridView<Grid> view(ebosSimulator_.vanguard().grid(), std::move(seeds[index]));
                domains_.emplace_back(index, std::move(partitions[index]), std::move(interior), std::move(view));
            }

            // Set up container for the local system matrices.
            domain_matrices_.resize(num_domains);

            // Set up container for the local linear solvers.
            for (int index = 0; index < num_domains; ++index) {
                // TODO: The ISTLSolverEbos constructor will make
                // parallel structures appropriate for the full grid
                // only. This must be addressed before going parallel.
                FlowLinearSolverParameters param;
                param.template init<TypeTag>();
                // Override solver type with umfpack if small domain.
                // Otherwise hardcode to ILU0
                if (domains_[index].cells.size() < 200) {
                    param.linsolver_ = "umfpack.json";
                } else {
                    param.linsolver_ = "ilu0";
                    param.linear_solver_reduction_ = 1e-2;
                }
                param.linear_solver_print_json_definition_ = false;
                domain_linsolvers_.emplace_back(ebosSimulator_, param);
            }

            assert(int(domains_.size()) == num_domains);
        }




        /// Called once before each time step.
        /// \param[in] timer                  simulation timer
        SimulatorReportSingle prepareStep(const SimulatorTimerInterface& timer)
        {
            SimulatorReportSingle report;
            Dune::Timer perfTimer;
            perfTimer.start();
            // update the solution variables in ebos
            if ( timer.lastStepFailed() ) {
                ebosSimulator_.model().updateFailed();
            } else {
                ebosSimulator_.model().advanceTimeLevel();
            }

            // Set the timestep size, episode index, and non-linear iteration index
            // for ebos explicitly. ebos needs to know the report step/episode index
            // because of timing dependent data despite the fact that flow uses its
            // own time stepper. (The length of the episode does not matter, though.)
            ebosSimulator_.setTime(timer.simulationTimeElapsed());
            ebosSimulator_.setTimeStepSize(timer.currentStepLength());
            ebosSimulator_.model().newtonMethod().setIterationIndex(0);

            ebosSimulator_.problem().beginTimeStep();

            unsigned numDof = ebosSimulator_.model().numGridDof();
            wasSwitched_.resize(numDof);
            std::fill(wasSwitched_.begin(), wasSwitched_.end(), false);

            if (param_.update_equations_scaling_) {
                std::cout << "equation scaling not supported yet" << std::endl;
                //updateEquationsScaling();
            }

            // Setup domain->well mapping.
            wellModel().setupDomains(domains_);

            report.pre_post_time += perfTimer.stop();

            return report;
        }


        /// Called once per nonlinear iteration.
        /// This model will perform a Newton-Raphson update, changing reservoir_state
        /// and well_state. It will also use the nonlinear_solver to do relaxation of
        /// updates if necessary.
        /// \param[in] iteration              should be 0 for the first call of a new timestep
        /// \param[in] timer                  simulation timer
        /// \param[in] nonlinear_solver       nonlinear solver used (for oscillation/relaxation control)
        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIteration(const int iteration,
                                                 const SimulatorTimerInterface& timer,
                                                 NonlinearSolverType& nonlinear_solver)
        {
            if (iteration == 0) {
                // For each iteration we store in a vector the norms of the residual of
                // the mass balance for each active phase, the well flux and the well equations.
                residual_norms_history_.clear();
                current_relaxation_ = 1.0;
                dx_old_ = 0.0;
                convergence_reports_.push_back({timer.reportStepNum(), timer.currentStepNum(), {}});
                convergence_reports_.back().report.reserve(11);
            }

            if (param_.nonlinear_solver_ == "aspin" || param_.nonlinear_solver_ == "nldd" || param_.nonlinear_solver_ == "purelocal") {
                return nonlinearIterationAspin(iteration, timer, nonlinear_solver);
            } else {
                return nonlinearIterationNewton(iteration, timer, nonlinear_solver);
            }
        }


        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIterationNewton(const int iteration,
                                                       const SimulatorTimerInterface& timer,
                                                       NonlinearSolverType& nonlinear_solver)
        {

            // -----------   Set up reports and timer   -----------
            SimulatorReportSingle report;
            failureReport_ = SimulatorReportSingle();
            Dune::Timer perfTimer;

            perfTimer.start();
            report.total_linearizations = 1;

            // -----------   Assemble   -----------
            try {
                report += assembleReservoir(timer, iteration);
                report.assemble_time += perfTimer.stop();
            }
            catch (...) {
                report.assemble_time += perfTimer.stop();
                failureReport_ += report;
                // todo (?): make the report an attribute of the class
                throw; // continue throwing the stick
            }

#if 0
            // DEBUG HACK
            // Check if domain-diagonal parts are identical to stored subdomain matrices.
            if (param_.nonlinear_solver_ == "nldd") {
                const Mat& main_matrix = ebosSimulator_.model().linearizer().jacobian().istlMatrix();
                for (const auto& domain : domains_) {
                    if (domain_matrices_[domain.index]) {
                        auto domain_diag_part = Details::extractMatrix(main_matrix, domain.cells);
                        if (Details::matrixEqual(domain_diag_part, *domain_matrices_[domain.index])) {
                            OpmLog::debug("  *** matrices IDENTICAL for domain " + std::to_string(domain.index));
                        } else {
                            OpmLog::debug("  *** matrices DIFFER for domain " + std::to_string(domain.index));
                        }
                    } else {
                        OpmLog::debug("  *** no local matrix for domain " + std::to_string(domain.index));
                    }
                }
            }
#endif

            // -----------   Check if converged   -----------
            std::vector<double> residual_norms;
            perfTimer.reset();
            perfTimer.start();
            // the step is not considered converged until at least minIter iterations is done
            {
                auto convrep = getConvergence(timer, iteration, residual_norms);
                // if (!convrep.converged()) {
                //     std::ostringstream os;
                //     os << "Convergence data for Newton iteration " << iteration << "\n"
                //        << convrep;
                //     OpmLog::debug(os.str());
                // }
                ////wellModel().logPrimaryVars();

                report.converged = convrep.converged()  && iteration > nonlinear_solver.minIter();;
                ConvergenceReport::Severity severity = convrep.severityOfWorstFailure();
                convergence_reports_.back().report.push_back(std::move(convrep));

                // Throw if any NaN or too large residual found.
                if (severity == ConvergenceReport::Severity::NotANumber) {
                    OPM_THROW(NumericalProblem, "NaN residual found!");
                } else if (severity == ConvergenceReport::Severity::TooLarge) {
                    OPM_THROW_NOLOG(NumericalProblem, "Too large residual found!");
                }
            }
            //report.convergence_check_time += perfTimer.stop();
            residual_norms_history_.push_back(residual_norms);
            // Dune::writeMatrixToMatlab(ebosSimulator_.model().linearizer().jacobian().istlMatrix(), "beforeWellModelLinearize.matlab");

            // -----------   If not converged, solve linear system   -----------
            if (!report.converged) {
                perfTimer.reset();
                perfTimer.start();
                report.total_newton_iterations = 1;

                // enable single precision for solvers when dt is smaller then 20 days
                //residual_.singlePrecision = (unit::convert::to(dt, unit::day) < 20.) ;

                // Compute the nonlinear update.
                unsigned nc = ebosSimulator_.model().numGridDof();
                BVector x(nc);

                // Solve the linear system.
                linear_solve_setup_time_ = 0.0;
                try {
                    // apply the Schur compliment of the well model to the reservoir linearized
                    // equations
                    // Note that linearize may throw for MSwells.
                    wellModel().linearize(ebosSimulator().model().linearizer().jacobian(),
                                          ebosSimulator().model().linearizer().residual());
                    //wellModel().logPrimaryVars();

                    solveJacobianSystem(x);
                    report.linear_solve_setup_time += linear_solve_setup_time_;
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();
                }
                catch (...) {
                    report.linear_solve_setup_time += linear_solve_setup_time_;
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();

                    failureReport_ += report;
                    throw; // re-throw up
                }

                perfTimer.reset();
                perfTimer.start();

                // handling well state update before oscillation treatment is a decision based
                // on observation to avoid some big performance degeneration under some circumstances.
                // there is no theorectical explanation which way is better for sure.
                wellModel().postSolve(x);
                //wellModel().logPrimaryVars();

                if (param_.use_update_stabilization_) {
                    // Stabilize the nonlinear update.
                    bool isOscillate = false;
                    bool isStagnate = false;
                    nonlinear_solver.detectOscillations(residual_norms_history_, residual_norms_history_.size() - 1, isOscillate, isStagnate);
                    if (isOscillate) {
                        current_relaxation_ -= nonlinear_solver.relaxIncrement();
                        current_relaxation_ = std::max(current_relaxation_, nonlinear_solver.relaxMax());
                        if (terminalOutputEnabled()) {
                            std::string msg = "    Oscillating behavior detected: Relaxation set to "
                                    + std::to_string(current_relaxation_);
                            OpmLog::info(msg);
                        }
                    }
                    nonlinear_solver.stabilizeNonlinearUpdate(x, dx_old_, current_relaxation_);
                }

                // Apply the update, with considering model-dependent limitations and
                // chopping of the update.
                updateSolution(x);

                report.update_time += perfTimer.stop();
            }

            return report;
        }




        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIterationAspin(const int iteration,
                                                      const SimulatorTimerInterface& timer,
                                                      NonlinearSolverType& nonlinear_solver)
        {
            // -----------   Set up reports and timer   -----------
            SimulatorReportSingle report;
            failureReport_ = SimulatorReportSingle();
            Dune::Timer perfTimer;

            //wellModel().logPrimaryVars();

            perfTimer.start();
            report.total_linearizations = 1;

            // -----------   Assemble   -----------
            try {
                report += assembleReservoir(timer, iteration);
                report.assemble_time += perfTimer.stop();
            }
            catch (...) {
                report.assemble_time += perfTimer.stop();
                failureReport_ += report;
                // todo (?): make the report an attribute of the class
                throw; // continue throwing the stick
            }
            //wellModel().logPrimaryVars();

            // -----------   Check if converged   -----------
            std::vector<double> residual_norms;
            perfTimer.reset();
            perfTimer.start();
            // the step is not considered converged until at least minIter iterations is done
            {
                auto convrep = getConvergence(timer, iteration, residual_norms);
                if (!convrep.converged()) {
                    // std::ostringstream os;
                    // os << "Convergence data for first global iteration:\n"
                    //    << convrep;
                    // OpmLog::debug(os.str());
                }
                report.converged = convrep.converged()  && iteration > nonlinear_solver.minIter();;
                ConvergenceReport::Severity severity = convrep.severityOfWorstFailure();
                convergence_reports_.back().report.push_back(std::move(convrep));

                // Throw if any NaN or too large residual found.
                if (severity == ConvergenceReport::Severity::NotANumber) {
                    OPM_THROW(NumericalProblem, "NaN residual found!");
                } else if (severity == ConvergenceReport::Severity::TooLarge) {
                    OPM_THROW_NOLOG(NumericalProblem, "Too large residual found!");
                }
            }
            //report.convergence_check_time += perfTimer.stop();
            residual_norms_history_.push_back(residual_norms);

            if (report.converged) {
                return report;
            }

            // -----------   If not converged, do an ASPIN or NLDD iteration   -----------

            // apply the Schur compliment of the well model to the reservoir linearized
            // equations
            // wellModel().linearize(ebosSimulator().model().linearizer().jacobian(),
            //                       ebosSimulator().model().linearizer().residual());

            // Take a copy of the FI residual.
            // auto fully_implicit_residual = ebosSimulator().model().linearizer().residual();
            // ... and the Jacobian.
            // auto fully_implicit_jacobian = ebosSimulator().model().linearizer().jacobian().istlMatrix();
            // ... and the solution state that generated it.
            auto& solution = ebosSimulator().model().solution(0);
            auto initial_solution = solution;
            // auto initial_wellstate = wellModel().wellState();
            auto locally_solved = initial_solution;

            // -----------   Decide on an ordering for the domains   -----------
            std::vector<int> domain_order(domains_.size());
            if (param_.local_solve_approach_ == "gauss-seidel") {
                if (true) {
                    // Use average pressures to order domains.
                    std::vector<std::pair<double, int>> avgpress_per_domain(domains_.size());
                    for (const auto& domain : domains_) {
                        double press_sum = 0.0;
                        for (const int c : domain.cells) {
                            press_sum += solution[c][Indices::pressureSwitchIdx];
                        }
                        const double avgpress = press_sum / domain.cells.size();
                        avgpress_per_domain[domain.index] = std::make_pair(avgpress, domain.index);
                    }
                    // Lexicographical sort by pressure, then index.
                    std::sort(avgpress_per_domain.begin(), avgpress_per_domain.end());
                    // Reverse since we want high-pressure regions solved first.
                    std::reverse(avgpress_per_domain.begin(), avgpress_per_domain.end());
                    for (size_t ii = 0; ii < domains_.size(); ++ii) {
                        domain_order[ii] = avgpress_per_domain[ii].second;
                    }
                } else {
                    // Use maximum residual to order domains.
                    const auto& residual = ebosSimulator_.model().linearizer().residual();
                    const int num_vars = residual[0].size();
                    std::vector<std::pair<double, int>> maxres_per_domain(domains_.size());
                    for (const auto& domain : domains_) {
                        double maxres = 0.0;
                        for (const int c : domain.cells) {
                            for (int ii = 0; ii < num_vars; ++ii) {
                                maxres = std::max(maxres, std::fabs(residual[c][ii]));
                            }
                        }
                        maxres_per_domain[domain.index] = std::make_pair(maxres, domain.index);
                    }
                    // Lexicographical sort by pressure, then index.
                    std::sort(maxres_per_domain.begin(), maxres_per_domain.end());
                    // Reverse since we want high-pressure regions solved first.
                    std::reverse(maxres_per_domain.begin(), maxres_per_domain.end());
                    for (size_t ii = 0; ii < domains_.size(); ++ii) {
                        domain_order[ii] = maxres_per_domain[ii].second;
                    }
                }
            } else {
                std::iota(domain_order.begin(), domain_order.end(), 0);
            }

            // -----------   Solve each domain separately   -----------
            std::vector<SimulatorReportSingle> domain_reports(domains_.size());
            for (const int domain_index : domain_order) {
                const auto& domain = domains_[domain_index];
                SimulatorReportSingle local_report;
                if (param_.local_solve_approach_ == "jacobi") {
                    auto initial_local_well_primary_vars = wellModel().getPrimaryVarsDomain(domain);
                    auto initial_local_solution = Details::extractVector(solution, domain.cells);
                    auto res = solveLocal(domain, timer, iteration);
                    local_report = res.first;
#if 0
                    auto local_solution = Details::extractVector(solution, domain.cells);
                    Details::setGlobal(local_solution, domain.cells, locally_solved);
                    Details::setGlobal(initial_local_solution, domain.cells, solution);
                   // ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain.view);
#else
                    if (/*!local_report.converged*/ false) {
                        // Try again with a less strict tolerance.
                        wellModel().setPrimaryVarsDomain(domain, initial_local_well_primary_vars);
                        Details::setGlobal(initial_local_solution, domain.cells, solution);
                        ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain.view);
                        // auto param_orig = param_;
                        // param_.local_tolerance_scaling_cnv_ *= 10.0;
                        // param_.local_tolerance_scaling_mb_ *= 10.0;
                        auto res = solveLocal(domain, timer, iteration);
                        local_report = res.first;
                        // param_ = param_orig;
                    }
                    if (local_report.converged) {
                        auto local_solution = Details::extractVector(solution, domain.cells);
                        Details::setGlobal(local_solution, domain.cells, locally_solved);
                        Details::setGlobal(initial_local_solution, domain.cells, solution);
                        ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain.view);
                    } else {
                        wellModel().setPrimaryVarsDomain(domain, initial_local_well_primary_vars);
                        Details::setGlobal(initial_local_solution, domain.cells, solution);
                        ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain.view);
                    }
#endif
                } else {
                    assert(param_.local_solve_approach_ == "gauss-seidel");
                    auto initial_local_well_primary_vars = wellModel().getPrimaryVarsDomain(domain);
                    auto initial_local_solution = Details::extractVector(solution, domain.cells);
                    auto res = solveLocal(domain, timer, iteration);
                    local_report = res.first;
                    if (/*!local_report.converged*/ false) {
                        // Try again with a less strict tolerance.
                        wellModel().setPrimaryVarsDomain(domain, initial_local_well_primary_vars);
                        Details::setGlobal(initial_local_solution, domain.cells, solution);
                        ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain.view);
                        // auto param_orig = param_;
                        // param_.local_tolerance_scaling_cnv_ = 50.0;
                        // param_.local_tolerance_scaling_mb_ = 50.0;
                        auto res = solveLocal(domain, timer, iteration);
                        local_report = res.first;
                        // param_ = param_orig;
                    }
                    if (!local_report.converged) {
                        // We look at the detailed convergence report to evaluate
                        // if we should accept the unconverged solution.
                        const auto& convrep = res.second;
                        // We do not accept a solution if the wells are unconverged.
                        if (!convrep.wellFailed()) {
                            // Calculare the sums of the mb and cnv failures.
                            double mb_sum = 0.0;
                            double cnv_sum = 0.0;
                            for (const auto& rf : convrep.reservoirFailures()) {
                                if (rf.type() == ConvergenceReport::ReservoirFailure::Type::MassBalance) {
                                    // mb_sum += rf.magnitude();
                                } else if (rf.type() == ConvergenceReport::ReservoirFailure::Type::Cnv) {
                                    // cnv_sum += rf.magnitude();
                                }
                            }
                            // If not too high, we overrule the convergence failure.
                            if (mb_sum < 1e-3 && cnv_sum < 1.0) {
                                local_report.converged = true;
                                OpmLog::debug("Accepting solution in unconverged domain " + std::to_string(domain.index));
                            }
                        }
                    }
                    if (local_report.converged) {
                        auto local_solution = Details::extractVector(solution, domain.cells);
                        Details::setGlobal(local_solution, domain.cells, locally_solved);
                    } else {
                        wellModel().setPrimaryVarsDomain(domain, initial_local_well_primary_vars);
                        Details::setGlobal(initial_local_solution, domain.cells, solution);
                        ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain.view);
                    }
                }
                // This should have updated the global matrix to be
                // dR_i/du_j evaluated at new local solutions for
                // i == j, at old solution for i != j.
                if (!local_report.converged) {
                    // TODO: more proper treatment, including in parallel.
                    OpmLog::debug("Convergence failure in domain " + std::to_string(domain.index));
                }
                domain_reports[domain.index] = local_report;
            }

#if 0
            // DEBUG write solution vectors.
            {
                const int nc = UgGridHelpers::numCells(grid_);
                BVector sol(nc);
                for (int c = 0; c < nc; ++c) {
                    sol[c] = solution[c];
                }
                BVector loc_sol(nc);
                for (int c = 0; c < nc; ++c) {
                    loc_sol[c] = locally_solved[c];
                }
                Dune::storeMatrixMarket(sol, param_.local_solve_approach_ + "-" + std::to_string(iteration) + ".mm");
                Dune::storeMatrixMarket(loc_sol, param_.local_solve_approach_ + "-local-" + std::to_string(iteration) + ".mm");
            }
#endif

            // Print summary of local solve convergence.
            {
                int num_converged = 0;
                SimulatorReportSingle rep;
                for (const auto& dr : domain_reports) {
                    if (dr.converged) {
                        ++num_converged;
                    }
                    rep += dr;
                }
                std::ostringstream os;
                os << fmt::format("Local solves finished. Converged for {}/{} domains.\n",
                                  num_converged, domain_reports.size());
                rep.reportFullyImplicit(os, nullptr);
                OpmLog::debug(os.str());
                local_reports_accumulated_ += rep;
            }

#define ASPIN_EXTRA_WELL_OUTPUT 1
#if ASPIN_EXTRA_WELL_OUTPUT
            // Extra debug output.
            //wellModel().logPrimaryVars();
            {
                std::ostringstream os;
                os << "   ** Well rates:";
                for (int w = 0; w < wellModel().wellState().numWells(); ++w) {
                    os << " | ";
                    const auto& wr = wellModel().wellState().wellRates(w);
                    for (const double r : wr) {
                        os << ' ' << r;
                    }
                }
                OpmLog::debug(os.str());
            }
#endif

            if (param_.nonlinear_solver_ == "purelocal") {
                if (param_.local_solve_approach_ == "jacobi") {
                    solution = locally_solved;
                    ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
                }
                return report;
            }
            else if (param_.nonlinear_solver_ == "nldd") {
                if (param_.local_solve_approach_ == "jacobi") {
                    solution = locally_solved;
                    ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
                }
                // Finish with a Newton step.
                // ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
                auto rep = nonlinearIterationNewton(iteration + 100, timer, nonlinear_solver);
                report += rep;
                if (rep.converged) {
                    report.converged = true;
                }
                return report;
            }
            report.total_newton_iterations = 1;

            /*
            // HACK to check FI convergence
            // Take a copy of the FI residual.
            auto res1 = ebosSimulator().model().linearizer().residual();
            // ... and the Jacobian.
            auto jac1 = ebosSimulator().model().linearizer().jacobian().istlMatrix();
            // ... and the solution state that generated it.
            auto sol1 = solution;
            auto ws1 = wellModel().wellState();

            // -----------   Assemble   -----------
            try {
                report += assembleReservoir(timer, iteration);
                report.assemble_time += perfTimer.stop();
            }
            catch (...) {
                report.assemble_time += perfTimer.stop();
                failureReport_ += report;
                // todo (?): make the report an attribute of the class
                throw; // continue throwing the stick
            }
            //wellModel().logPrimaryVars();

            // -----------   Check if converged   -----------
            perfTimer.reset();
            perfTimer.start();
            // the step is not considered converged until at least minIter iterations is done
            {
                auto convrep = getConvergence(timer, iteration, residual_norms);
                report.converged = convrep.converged()  && iteration > nonlinear_solver.minIter();;
                ConvergenceReport::Severity severity = convrep.severityOfWorstFailure();
                convergence_reports_.back().report.push_back(std::move(convrep));

                // Throw if any NaN or too large residual found.
                if (severity == ConvergenceReport::Severity::NotANumber) {
                    OPM_THROW(NumericalProblem, "NaN residual found!");
                } else if (severity == ConvergenceReport::Severity::TooLarge) {
                    OPM_THROW(NumericalProblem, "Too large residual found!");
                }
            }
            report.update_time += perfTimer.stop();
            residual_norms_history_.push_back(residual_norms);

            if (report.converged) {
                return report;
            }
            // HACK reinstate everything
            ebosSimulator().model().linearizer().residual() = res1;
            ebosSimulator().model().linearizer().jacobian().istlMatrix() = jac1;
            solution = sol1;
            wellModel().wellState() = ws1;
            */

            // -----------   Compute ASPIN residual, check convergence   -----------
            const int nc = grid_.size(0); // Includes both owned cells and overlap.
            BVector aspin_residual(nc);
            const int num_vars = aspin_residual[0].size();
            std::vector<double> max_sol(num_vars, 0.0);
            std::vector<double> max_diff(num_vars, 0.0);
            for (int c = 0; c < nc; ++c) {
                aspin_residual[c] = initial_solution[c] - locally_solved[c];
                for (int var = 0; var < num_vars; ++var) {
                    max_sol[var] = std::max(max_sol[var], std::fabs(locally_solved[c][var]));
                    max_diff[var] = std::max(max_diff[var], std::fabs(aspin_residual[c][var]));
                }
            }
            double max_scaled_diff = 0.0;
            for (int var = 0; var < num_vars; ++var) {
                max_scaled_diff = std::max(max_scaled_diff, max_diff[var]/max_sol[var]);
            }
            OpmLog::debug("  Scaled ASPIN residual: " + std::to_string(max_scaled_diff));
            if (max_scaled_diff < param_.outer_aspin_tolerance_) {
                report.converged = true;
                solution = locally_solved;
                return report;
            }

            // -----------   Solve global linear system   -----------
            // The matrix is kept as is for now.
            // TODO: test with proper ASPIN (in what we have now, columns are updated outside the diagonal
            // domain-domain interaction blocks in a Gauss-Seidel-like manner), and other variations.
            // NOTE: the above should no longer be true, the matrix is ASPIN now (should be).
            //
#if 0
            // Compute rhs.
            // The modified residual is D * aspin_residual in each domain.
            auto& ebosResid = ebosSimulator_.model().linearizer().residual();
            for (const auto& domain : domains_) {
                auto aspin_res_local = Details::extractVector(aspin_residual, domain.cells);
                BVector aspin_res_mod(aspin_res_local.size());
                aspin_res_mod = 0;
                if (domain_reports[domain.index].total_newton_iterations > 0) {
                    (*domain_matrices_[domain.index]).mv(aspin_res_local, aspin_res_mod);
                }
                Details::setGlobal(aspin_res_mod, domain.cells, ebosResid);
            }
            // Add J_aspin * aspin_residual to the modified residual, so we solve for
            // the difference from the local solution rather than from the very start
            // of the iteration.
            {
                const auto& ebosJac = ebosSimulator_.model().linearizer().jacobian().istlMatrix();
                BVector yy(nc);
                ebosJac.mv(aspin_residual, yy);
                ebosResid -= yy;
            }
#else
            // Compute rhs.
            {
                auto Jcopy = ebosSimulator_.model().linearizer().jacobian().istlMatrix();
                Dune::storeMatrixMarket(Jcopy, "aspinjac.mm");
                // Eliminate the domain diagonal blocks from Jcopy. Jcopy is then (J - D)
                for (const auto& domain : domains_) {
                    for (const int cell : domain.cells) {
                        auto& row = Jcopy[cell];
                        for (auto iter = row.begin(); iter != row.end(); ++iter) {
                            int col = iter.index();
                            auto lb = std::find(domain.cells.begin(), domain.cells.end(), col);
                            if (lb != domain.cells.end()) {
                                //std::cerr << "row = " << cell << "   col = " << col <<'\n';
                                *iter = 0.0;
                            }
                        }
                    }
                }
                Dune::storeMatrixMarket(Jcopy, "aspinjacoffdiag.mm");
                // Set the residual to -Jcopy * aspin_residual, which will then be (D - J)*(x^k - x^{k+1})
                auto& ebosResid = ebosSimulator_.model().linearizer().residual();
                ebosResid = 0.0;
                Jcopy.mmv(aspin_residual, ebosResid);
            }
#endif

            // Check convergence before solving.
            {
                auto convrep = getConvergence(timer, -iteration, residual_norms);
                // report.converged = convrep.converged()  && iteration > nonlinear_solver.minIter();;
                // ConvergenceReport::Severity severity = convrep.severityOfWorstFailure();
                // convergence_reports_.back().report.push_back(std::move(convrep));

                // // Throw if any NaN or too large residual found.
                // if (severity == ConvergenceReport::Severity::NotANumber) {
                //     OPM_THROW(NumericalProblem, "NaN residual found!");
                // } else if (severity == ConvergenceReport::Severity::TooLarge) {
                //     OPM_THROW_NOLOG(NumericalProblem, "Too large residual found!");
                // }
            }
            // report.convergence_check_time += perfTimer.stop();
            // residual_norms_history_.push_back(residual_norms);

            // if (report.converged) {
            //     return report;
            // }


            // Solve the linear system.
            linear_solve_setup_time_ = 0.0;
            BVector x(nc);
            try {
                solveJacobianSystem(x);
                report.linear_solve_setup_time += linear_solve_setup_time_;
                report.linear_solve_time += perfTimer.stop();
                report.total_linear_iterations += linearIterationsLastSolve();
            }
            catch (...) {
                report.linear_solve_setup_time += linear_solve_setup_time_;
                report.linear_solve_time += perfTimer.stop();
                report.total_linear_iterations += linearIterationsLastSolve();

                failureReport_ += report;
                throw; // re-throw up
            }

            perfTimer.reset();
            perfTimer.start();

            // -----------   Update solution   -----------

            // Compute the nonlinear update.

            // apply the Schur compliment of the well model to the reservoir linearized
            // equations
            // wellModel().linearize(ebosSimulator().model().linearizer().jacobian(),
            //                       ebosSimulator().model().linearizer().residual());


            // handling well state update before oscillation treatment is a decision based
            // on observation to avoid some big performance degeneration under some circumstances.
            // there is no theorectical explanation which way is better for sure.
            // wellModel().wellState() = initial_wellstate;
            // ebosSimulator_.model().newtonMethod().setIterationIndex(iteration);
            // wellModel().beginIteration();

            // Before we can do postSolve() we must take into account that the current
            // state in the wells are from after the local solves, not from the start
            // of the iteration.
            // Instead of passing x (where state_{i+1} = state_i + x, i.e. x is the
            // difference from one global iteration to the next) we should pass
            // y, with y = state_{i+1} - state_after_local_solve.
            // BVector y(nc);
            // for (int c = 0; c < nc; ++c) {
            //     y[c] = solution[c] - locally_solved[c];
            // }
            // wellModel().postSolve(y);

            if (param_.use_update_stabilization_) {
                // Stabilize the nonlinear update.
                bool isOscillate = false;
                bool isStagnate = false;
                nonlinear_solver.detectOscillations(residual_norms_history_, iteration, isOscillate, isStagnate);
                if (isOscillate) {
                    current_relaxation_ -= nonlinear_solver.relaxIncrement();
                    current_relaxation_ = std::max(current_relaxation_, nonlinear_solver.relaxMax());
                    if (terminalOutputEnabled()) {
                        std::string msg = "    Oscillating behavior detected: Relaxation set to "
                            + std::to_string(current_relaxation_);
                        OpmLog::info(msg);
                    }
                }
                nonlinear_solver.stabilizeNonlinearUpdate(x, dx_old_, current_relaxation_);
            }

            // Apply the update, with considering model-dependent limitations and
            // chopping of the update. This update must be based on the initial solution,
            // not the solution after the local solve.
            // ebosSimulator().model().solution(0) = initial_solution;
            updateSolution(x);




            // Before we can do postSolve() we must take into account that the current
            // state in the wells are from after the local solves, not from the start
            // of the iteration.
            // Instead of passing x (where state_{i+1} = state_i - x, i.e. x is the
            // difference from one global iteration to the next) we should pass
            // y, with y = state_after_local_solve - state_{i+1}.
            // BVector y(nc);
            // for (int c = 0; c < nc; ++c) {
            //     y[c] = locally_solved[c] - solution[c];
            // }
            // wellModel().postSolve(y);
            wellModel().postSolve(x);




            report.update_time += perfTimer.stop();

            //wellModel().logPrimaryVars();

            // Extra debug output.
            {
                std::ostringstream os;
                os << "   ** Well rates:";
                for (int w = 0; w < wellModel().wellState().numWells(); ++w) {
                    os << " | ";
                    const auto& wr = wellModel().wellState().wellRates(w);
                    for (const double r : wr) {
                        os << ' ' << r;
                    }
                }
                OpmLog::debug(os.str());
            }


            return report;
        }




        std::pair<SimulatorReportSingle, ConvergenceReport>
        solveLocal(const Domain& domain,
                   const SimulatorTimerInterface& timer,
                   const int global_iteration,
                   const bool initial_assembly_required = false)
        {
            SimulatorReportSingle report;
            Dune::Timer solveTimer;
            solveTimer.start();
            Dune::Timer detailTimer;

            // TODO: the following assertion will fail, does it indicate a programming error?
            // assert(ebosSimulator_.model().newtonMethod().numIterations() == 0);
            ebosSimulator_.model().newtonMethod().setIterationIndex(0);
            //ebosSimulator_.problem().beginIteration();

            // When called, if assembly has already been performed
            // with the initial values, we only need to check
            // for local convergence. Otherwise, we must do a local
            // assembly.
            int iter = 0;
            if (initial_assembly_required) {
                OpmLog::debug("solveLocal() initial assembly START");
                //wellModel().logPrimaryVars();
                detailTimer.start();
                ebosSimulator_.model().newtonMethod().setIterationIndex(iter);
                // TODO: we should have a beginIterationLocal function()
                // only handling the well model for now
                // ebosSimulator_.problem().beginIteration();
                ebosSimulator_.problem().wellModel().assembleLocal(ebosSimulator_.model().newtonMethod().numIterations(),
                                                                   ebosSimulator_.timeStepSize(),
                                                                   domain);
                // Assemble reservoir locally.
                report += assembleReservoirLocal(domain, iter);
                report.assemble_time += detailTimer.stop();
                OpmLog::debug("solveLocal() initial assembly END");
                //wellModel().logPrimaryVars();
            }
            detailTimer.reset();
            detailTimer.start();
            std::vector<double> resnorms;
            auto convreport = getLocalConvergence(domain, timer, 0, resnorms);
            if (convreport.converged()) {
                // TODO: set more info, timing etc.
                report.converged = true;
                // Exit, but not if we are on the first global iteration
                // (if so, we require at least one iteration).
                //if (global_iteration > 0) {
                    return { report, convreport };
                //}
            } else {
                // std::ostringstream os;
                // os << "Convergence data for first local iteration:\n"
                //    << convreport;
                // OpmLog::debug(os.str());
            }
            //report.convergence_check_time += detailTimer.stop();

            // We have already assembled for the first iteration,
            // but not done the Schur complement for the wells yet.
            detailTimer.reset();
            detailTimer.start();
            wellModel().linearizeDomain(domain,
                                        ebosSimulator().model().linearizer().jacobian(),
                                        ebosSimulator().model().linearizer().residual());
            const double tt1 = detailTimer.stop();
            report.assemble_time += tt1;
            report.assemble_time_well += tt1;
            //wellModel().logPrimaryVars();

            // Local Newton loop.
            const int max_iter = param_.max_local_solve_iterations_;
            bool converged = false;
            do {
                // Solve local linear system.
                // Note that x has full size, we expect it to be nonzero only for in-domain cells.
                const int nc = grid_.size(0);
                BVector x(nc);
                detailTimer.reset();
                detailTimer.start();
                solveLocalJacobianSystem(domain, x);
                wellModel().postSolveLocal(x, domain);
                report.linear_solve_time += detailTimer.stop();
                report.linear_solve_setup_time += linear_solve_setup_time_;
                report.total_linear_iterations = linearIterationsLastSolve();
                //wellModel().logPrimaryVars();

                // Update local solution. // TODO: x is still full size, should we optimize it?
                detailTimer.reset();
                detailTimer.start();
                updateDomainSolution(domain, x);
                report.update_time += detailTimer.stop();
                // Note: endIteration does not do anything
                // ebosSimulator_.problem().endIteration();

                // Assemble well and reservoir.
                detailTimer.reset();
                detailTimer.start();
                ++iter;
                ebosSimulator_.model().newtonMethod().setIterationIndex(iter);
                // TODO: we should have a beginIterationLocal function()
                // only handling the well model for now
                // ebosSimulator_.problem().beginIteration();
                // Assemble reservoir locally.
                ebosSimulator_.problem().wellModel().assembleLocal(ebosSimulator_.model().newtonMethod().numIterations(),
                                                                   ebosSimulator_.timeStepSize(),
                                                                   domain);
                report += assembleReservoirLocal(domain, iter);
                report.assemble_time += detailTimer.stop();


                // Check for local convergence.
                detailTimer.reset();
                detailTimer.start();
                convreport = getLocalConvergence(domain, timer, iter, resnorms);
                // report.convergence_check_time += detailTimer.stop();
                //wellModel().logPrimaryVars();
                // Dune::writeMatrixToMatlab(ebosSimulator_.model().linearizer().jacobian().istlMatrix(), "beforeWellModelLinearize.matlab");

                // apply the Schur compliment of the well model to the reservoir linearized
                // equations
                detailTimer.reset();
                detailTimer.start();
                wellModel().linearizeDomain(domain,
                                            ebosSimulator().model().linearizer().jacobian(),
                                            ebosSimulator().model().linearizer().residual());
                const double tt2 = detailTimer.stop();
                report.assemble_time += tt2;
                report.assemble_time_well += tt2;
                //wellModel().logPrimaryVars();

            } while (!convreport.converged() && iter <= max_iter);

            ebosSimulator_.problem().endIteration();

            if (!convreport.converged()) {
                // std::ostringstream os;
                // os << "Convergence data for unconverged local solution:\n"
                //    << convreport;
                // OpmLog::debug(os.str());
            }

            report.converged = convreport.converged();
            report.total_newton_iterations = iter;
            report.total_linearizations = iter;
            report.total_time = solveTimer.stop();
            // TODO: set more info, timing etc.
            return { report, convreport };
        }




        void printIf(int c, double x, double y, double eps, std::string type) {
            if (std::abs(x-y) > eps) {
                std::cout << type << " " <<c << ": "<<x << " " << y << std::endl;
            }
        }




        /// Called once after each time step.
        /// In this class, this function does nothing.
        /// \param[in] timer                  simulation timer
        SimulatorReportSingle afterStep(const SimulatorTimerInterface&)
        {
            SimulatorReportSingle report;
            Dune::Timer perfTimer;
            perfTimer.start();
            ebosSimulator_.problem().endTimeStep();
            report.pre_post_time += perfTimer.stop();
            // report.completed_timesteps = 1;
            return report;
        }

        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        SimulatorReportSingle assembleReservoir(const SimulatorTimerInterface& /* timer */,
                                                const int iterationIdx)
        {
            // TODO: Undo this temporary testing change.
            // if (param_.enable_aspin_) {
            //     ebosSimulator_.model().newtonMethod().setIterationIndex(iterationIdx);
            //     ebosSimulator_.problem().beginIteration();
            //     ebosSimulator_.model().linearizer().resetSystem();
            //     for (const auto& domain : domains_) {
            //         assembleReservoirLocal(domain, iterationIdx);
            //     }
            //     ebosSimulator_.problem().endIteration();
            // } else {
                ebosSimulator_.model().newtonMethod().setIterationIndex(iterationIdx);
                ebosSimulator_.problem().beginIteration();
                ebosSimulator_.model().linearizer().linearizeDomain();
                ebosSimulator_.problem().endIteration();
            // }
            return wellModel().lastReport();
        }

        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        SimulatorReportSingle assembleReservoirLocal(const Domain& domain,
                                                     const int iterationIdx)
        {
            // -------- Mass balance equations --------
            // ebosSimulator_.model().newtonMethod().setIterationIndex(iterationIdx);
            // ebosSimulator_.problem().beginIteration();
            // Need to set residual and jacobian to zero in the domain.
            ebosSimulator_.model().linearizer().resetSystem(domain);
            // Call the domain-dependent linearization.
            ebosSimulator_.model().linearizer().linearizeDomain(domain);
            // ebosSimulator_.problem().endIteration();

            return wellModel().lastReport();
        }

        // compute the "relative" change of the solution between time steps
        double relativeChange() const
        {
            Scalar resultDelta = 0.0;
            Scalar resultDenom = 0.0;

            const auto& elemMapper = ebosSimulator_.model().elementMapper();
            const auto& gridView = ebosSimulator_.gridView();
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                unsigned globalElemIdx = elemMapper.index(elem);
                const auto& priVarsNew = ebosSimulator_.model().solution(/*timeIdx=*/0)[globalElemIdx];

                Scalar pressureNew;
                pressureNew = priVarsNew[Indices::pressureSwitchIdx];

                Scalar saturationsNew[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationNew = 1.0;
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
                    FluidSystem::numActivePhases() > 1 &&
                    priVarsNew.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                    saturationsNew[FluidSystem::waterPhaseIdx] = priVarsNew[Indices::waterSwitchIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::waterPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                    FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                    priVarsNew.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg) {
                    assert(Indices::compositionSwitchIdx >= 0 );
                    saturationsNew[FluidSystem::gasPhaseIdx] = priVarsNew[Indices::compositionSwitchIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::gasPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    saturationsNew[FluidSystem::oilPhaseIdx] = oilSaturationNew;
                }

                const auto& priVarsOld = ebosSimulator_.model().solution(/*timeIdx=*/1)[globalElemIdx];

                Scalar pressureOld;
                pressureOld = priVarsOld[Indices::pressureSwitchIdx];

                Scalar saturationsOld[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationOld = 1.0;

                // NB fix me! adding pressures changes to satutation changes does not make sense
                Scalar tmp = pressureNew - pressureOld;
                resultDelta += tmp*tmp;
                resultDenom += pressureNew*pressureNew;

                if (FluidSystem::numActivePhases() > 1) {
                    if (priVarsOld.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                        saturationsOld[FluidSystem::waterPhaseIdx] = priVarsOld[Indices::waterSwitchIdx];
                        oilSaturationOld -= saturationsOld[FluidSystem::waterPhaseIdx];
                    }

                    if (priVarsOld.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg)
                    {
                        assert(Indices::compositionSwitchIdx >= 0 );
                        saturationsOld[FluidSystem::gasPhaseIdx] = priVarsOld[Indices::compositionSwitchIdx];
                        oilSaturationOld -= saturationsOld[FluidSystem::gasPhaseIdx];
                    }

                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        saturationsOld[FluidSystem::oilPhaseIdx] = oilSaturationOld;
                    }
                    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++ phaseIdx) {
                        Scalar tmpSat = saturationsNew[phaseIdx] - saturationsOld[phaseIdx];
                        resultDelta += tmpSat*tmpSat;
                        resultDenom += saturationsNew[phaseIdx]*saturationsNew[phaseIdx];
                        assert(std::isfinite(resultDelta));
                        assert(std::isfinite(resultDenom));
                    }
                }
            }

            resultDelta = gridView.comm().sum(resultDelta);
            resultDenom = gridView.comm().sum(resultDenom);

            if (resultDenom > 0.0)
                return resultDelta/resultDenom;
            return 0.0;
        }


        /// Number of linear iterations used in last call to solveJacobianSystem().
        int linearIterationsLastSolve() const
        {
            return ebosSimulator_.model().newtonMethod().linearSolver().iterations ();
        }


        void solveLocalJacobianSystem(const Domain& domain, BVector& global_x)
        {
            Dune::Timer perfTimer;
            perfTimer.start();

            const Mat& main_matrix = ebosSimulator_.model().linearizer().jacobian().istlMatrix();
            if (domain_matrices_[domain.index]) {
                Details::copySubMatrix(main_matrix, domain.cells, *domain_matrices_[domain.index]);
            } else {
                domain_matrices_[domain.index] = std::make_unique<Mat>(Details::extractMatrix(main_matrix, domain.cells));
            }
            auto& jac = *domain_matrices_[domain.index];
            auto res = Details::extractVector(ebosSimulator_.model().linearizer().residual(), domain.cells);
            auto x = res;

            // set initial guess
            global_x = 0.0;
            x = 0.0;

            auto& linsolver = domain_linsolvers_[domain.index];

            linsolver.prepare(jac, res);
            linear_solve_setup_time_ = perfTimer.stop();
            linsolver.setResidual(res);
            //linsolver.setMatrix(jac);
            linsolver.solve(x);

            Details::setGlobal(x, domain.cells, global_x);
        }

        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        void solveJacobianSystem(BVector& x)
        {

            auto& ebosJac = ebosSimulator_.model().linearizer().jacobian().istlMatrix();
            auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            // set initial guess
            x = 0.0;

            auto& ebosSolver = ebosSimulator_.model().newtonMethod().linearSolver();
            Dune::Timer perfTimer;
            perfTimer.start();
            ebosSolver.prepare(ebosJac, ebosResid);
            linear_solve_setup_time_ = perfTimer.stop();
            ebosSolver.setResidual(ebosResid);
            // actually, the error needs to be calculated after setResidual in order to
            // account for parallelization properly. since the residual of ECFV
            // discretizations does not need to be synchronized across processes to be
            // consistent, this is not relevant for OPM-flow...
            //ebosSolver.setMatrix(ebosJac);
            ebosSolver.solve(x);
       }



        /// Apply an update to the primary variables.
        void updateDomainSolution(const Domain& domain, const BVector& dx)
        {
            auto& ebosNewtonMethod = ebosSimulator_.model().newtonMethod();
            SolutionVector& solution = ebosSimulator_.model().solution(/*timeIdx=*/0);

            ebosNewtonMethod.update_(/*nextSolution=*/solution,
                                     /*curSolution=*/solution,
                                     /*update=*/dx,
                                     /*resid=*/dx,
                                     domain.cells); // the update routines of the black
                                                    // oil model do not care about the
                                                    // residual

            // if the solution is updated, the intensive quantities need to be recalculated
            ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain.view);
        }




        /// Apply an update to the primary variables.
        void updateSolution(const BVector& dx)
        {
            OPM_TIMEBLOCK(updateSolution);
            auto& ebosNewtonMethod = ebosSimulator_.model().newtonMethod();
            SolutionVector& solution = ebosSimulator_.model().solution(/*timeIdx=*/0);

            ebosNewtonMethod.update_(/*nextSolution=*/solution,
                                     /*curSolution=*/solution,
                                     /*update=*/dx,
                                     /*resid=*/dx); // the update routines of the black
                                                    // oil model do not care about the
                                                    // residual

            // if the solution is updated, the intensive quantities need to be recalculated
            {
                OPM_TIMEBLOCK(invalidateAndUpdateIntensiveQuantities);
                ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
                ebosSimulator_.problem().eclWriter()->mutableEclOutputModule().invalidateLocalData();
            }
        }

        /// Return true if output to cout is wanted.
        bool terminalOutputEnabled() const
        {
            return terminal_output_;
        }

        template <class CollectiveCommunication>
        std::tuple<double,double> convergenceReduction(const CollectiveCommunication& comm,
                                                       const double pvSumLocal,
                                                       const double numAquiferPvSumLocal,
                                                       std::vector< Scalar >& R_sum,
                                                       std::vector< Scalar >& maxCoeff,
                                                       std::vector< Scalar >& B_avg)
        {
            OPM_TIMEBLOCK(convergenceReduction);
            // Compute total pore volume (use only owned entries)
            double pvSum = pvSumLocal;
            double numAquiferPvSum = numAquiferPvSumLocal;

            if( comm.size() > 1 )
            {
                // global reduction
                std::vector< Scalar > sumBuffer;
                std::vector< Scalar > maxBuffer;
                const int numComp = B_avg.size();
                sumBuffer.reserve( 2*numComp + 2 ); // +2 for (numAquifer)pvSum
                maxBuffer.reserve( numComp );
                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    sumBuffer.push_back( B_avg[ compIdx ] );
                    sumBuffer.push_back( R_sum[ compIdx ] );
                    maxBuffer.push_back( maxCoeff[ compIdx ] );
                }

                // Compute total pore volume
                sumBuffer.push_back( pvSum );
                sumBuffer.push_back( numAquiferPvSum );

                // compute global sum
                comm.sum( sumBuffer.data(), sumBuffer.size() );

                // compute global max
                comm.max( maxBuffer.data(), maxBuffer.size() );

                // restore values to local variables
                for( int compIdx = 0, buffIdx = 0; compIdx < numComp; ++compIdx, ++buffIdx )
                {
                    B_avg[ compIdx ]    = sumBuffer[ buffIdx ];
                    ++buffIdx;

                    R_sum[ compIdx ]       = sumBuffer[ buffIdx ];
                }

                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    maxCoeff[ compIdx ] = maxBuffer[ compIdx ];
                }

                // restore global pore volume
                pvSum = sumBuffer[sumBuffer.size()-2];
                numAquiferPvSum = sumBuffer.back();
            }

            // return global pore volume
            return {pvSum, numAquiferPvSum};
        }

        /// \brief Get reservoir quantities on this process needed for convergence calculations.
        /// \return A pair of the local pore volume of interior cells and the pore volumes
        ///         of the cells associated with a numerical aquifer.
        std::tuple<double,double> localConvergenceData(std::vector<Scalar>& R_sum,
                                    std::vector<Scalar>& maxCoeff,
                                    std::vector<Scalar>& B_avg,
                                    std::vector<int>& maxCoeffCell)
        {
            OPM_TIMEBLOCK(localConvergenceData);
            double pvSumLocal = 0.0;
            double numAquiferPvSumLocal = 0.0;
            const auto& ebosModel = ebosSimulator_.model();
            const auto& ebosProblem = ebosSimulator_.problem();

            const auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            ElementContext elemCtx(ebosSimulator_);
            const auto& gridView = ebosSimulator().gridView();
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                elemCtx.updatePrimaryStencil(elem);
                // elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                // const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
                const auto& fs = intQuants.fluidState();

                const double pvValue = ebosProblem.referencePorosity(cell_idx, /*timeIdx=*/0) * ebosModel.dofTotalVolume( cell_idx );
                pvSumLocal += pvValue;

                if (isNumericalAquiferCell(gridView.grid(), elem))
                {
                    numAquiferPvSumLocal += pvValue;
                }

                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }

                    const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));

                    B_avg[ compIdx ] += 1.0 / fs.invB(phaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][compIdx];

                    R_sum[ compIdx ] += R2;
                    const double Rval = std::abs( R2 ) / pvValue;
                    if (Rval > maxCoeff[ compIdx ]) {
                        maxCoeff[ compIdx ] = Rval;
                        maxCoeffCell[ compIdx ] = cell_idx;
                    }
                }

                if constexpr (has_solvent_) {
                    B_avg[ contiSolventEqIdx ] += 1.0 / intQuants.solventInverseFormationVolumeFactor().value();
                    const auto R2 = ebosResid[cell_idx][contiSolventEqIdx];
                    R_sum[ contiSolventEqIdx ] += R2;
                    maxCoeff[ contiSolventEqIdx ] = std::max( maxCoeff[ contiSolventEqIdx ], std::abs( R2 ) / pvValue );
                }
                if constexpr (has_extbo_) {
                    B_avg[ contiZfracEqIdx ] += 1.0 / fs.invB(FluidSystem::gasPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiZfracEqIdx];
                    R_sum[ contiZfracEqIdx ] += R2;
                    maxCoeff[ contiZfracEqIdx ] = std::max( maxCoeff[ contiZfracEqIdx ], std::abs( R2 ) / pvValue );
                }
                if constexpr (has_polymer_) {
                    B_avg[ contiPolymerEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiPolymerEqIdx];
                    R_sum[ contiPolymerEqIdx ] += R2;
                    maxCoeff[ contiPolymerEqIdx ] = std::max( maxCoeff[ contiPolymerEqIdx ], std::abs( R2 ) / pvValue );
                }
                if constexpr (has_foam_) {
                    B_avg[ contiFoamEqIdx ] += 1.0 / fs.invB(FluidSystem::gasPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiFoamEqIdx];
                    R_sum[ contiFoamEqIdx ] += R2;
                    maxCoeff[ contiFoamEqIdx ] = std::max( maxCoeff[ contiFoamEqIdx ], std::abs( R2 ) / pvValue );
                }
                if constexpr (has_brine_) {
                    B_avg[ contiBrineEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiBrineEqIdx];
                    R_sum[ contiBrineEqIdx ] += R2;
                    maxCoeff[ contiBrineEqIdx ] = std::max( maxCoeff[ contiBrineEqIdx ], std::abs( R2 ) / pvValue );
                }

                if constexpr (has_polymermw_) {
                    static_assert(has_polymer_);

                    B_avg[contiPolymerMWEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    // the residual of the polymer molecular equation is scaled down by a 100, since molecular weight
                    // can be much bigger than 1, and this equation shares the same tolerance with other mass balance equations
                    // TODO: there should be a more general way to determine the scaling-down coefficient
                    const auto R2 = ebosResid[cell_idx][contiPolymerMWEqIdx] / 100.;
                    R_sum[contiPolymerMWEqIdx] += R2;
                    maxCoeff[contiPolymerMWEqIdx] = std::max( maxCoeff[contiPolymerMWEqIdx], std::abs( R2 ) / pvValue );
                }

                if constexpr (has_energy_) {
                    B_avg[ contiEnergyEqIdx ] += 1.0;
                    const auto R2 = ebosResid[cell_idx][contiEnergyEqIdx];
                    R_sum[ contiEnergyEqIdx ] += R2;
                    maxCoeff[ contiEnergyEqIdx ] = std::max( maxCoeff[ contiEnergyEqIdx ], std::abs( R2 ) / pvValue );
                }

            }

            OPM_END_PARALLEL_TRY_CATCH("BlackoilModelEbos::localConvergenceData() failed: ", grid_.comm());

            // compute local average in terms of global number of elements
            const int bSize = B_avg.size();
            for ( int i = 0; i<bSize; ++i )
            {
                B_avg[ i ] /= Scalar( global_nc_ );
            }

            return {pvSumLocal, numAquiferPvSumLocal};
        }


        // Get reservoir quantities on this process needed for convergence calculations.
        double localDomainConvergenceData(const Domain& domain,
                                          std::vector<Scalar>& R_sum,
                                          std::vector<Scalar>& maxCoeff,
                                          std::vector<Scalar>& B_avg,
                                          std::vector<int>& maxCoeffCell)
        {
            double pvSumLocal = 0.0;
            const auto& ebosModel = ebosSimulator_.model();
            const auto& ebosProblem = ebosSimulator_.problem();

            const auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            ElementContext elemCtx(ebosSimulator_);
            const auto& gridView = domain.view;
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            for (auto elemIt = gridView.template begin</*codim=*/0>();
                 elemIt != elemEndIt;
                 ++elemIt)
            {
                if (elemIt->partitionType() != Dune::InteriorEntity) {
                    continue;
                }
                const auto& elem = *elemIt;
                elemCtx.updatePrimaryStencil(elem);
                // elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);

                // const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                auto* intQuantsPtr = ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0);
                // const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
                if (intQuantsPtr == nullptr) {
                    OpmLog::debug("** no cached intensive quantities for cell " + std::to_string(cell_idx));
                    elemCtx.updatePrimaryIntensiveQuantities(0);
                    intQuantsPtr = ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0);
                    assert(intQuantsPtr != nullptr);
                }
                const auto& intQuants = *intQuantsPtr;

                const auto& fs = intQuants.fluidState();

                const double pvValue = ebosProblem.referencePorosity(cell_idx, /*timeIdx=*/0) * ebosModel.dofTotalVolume( cell_idx );
                pvSumLocal += pvValue;

                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }

                    const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));

                    B_avg[ compIdx ] += 1.0 / fs.invB(phaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][compIdx];

                    R_sum[ compIdx ] += R2;
                    const double Rval = std::abs( R2 ) / pvValue;
                    if (Rval > maxCoeff[ compIdx ]) {
                        maxCoeff[ compIdx ] = Rval;
                        maxCoeffCell[ compIdx ] = cell_idx;
                    }
                }

                if constexpr (has_solvent_) {
                    B_avg[ contiSolventEqIdx ] += 1.0 / intQuants.solventInverseFormationVolumeFactor().value();
                    const auto R2 = ebosResid[cell_idx][contiSolventEqIdx];
                    R_sum[ contiSolventEqIdx ] += R2;
                    maxCoeff[ contiSolventEqIdx ] = std::max( maxCoeff[ contiSolventEqIdx ], std::abs( R2 ) / pvValue );
                }
                if constexpr (has_extbo_) {
                    B_avg[ contiZfracEqIdx ] += 1.0 / fs.invB(FluidSystem::gasPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiZfracEqIdx];
                    R_sum[ contiZfracEqIdx ] += R2;
                    maxCoeff[ contiZfracEqIdx ] = std::max( maxCoeff[ contiZfracEqIdx ], std::abs( R2 ) / pvValue );
                }
                if constexpr (has_polymer_) {
                    B_avg[ contiPolymerEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiPolymerEqIdx];
                    R_sum[ contiPolymerEqIdx ] += R2;
                    maxCoeff[ contiPolymerEqIdx ] = std::max( maxCoeff[ contiPolymerEqIdx ], std::abs( R2 ) / pvValue );
                }
                if constexpr (has_foam_) {
                    B_avg[ contiFoamEqIdx ] += 1.0 / fs.invB(FluidSystem::gasPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiFoamEqIdx];
                    R_sum[ contiFoamEqIdx ] += R2;
                    maxCoeff[ contiFoamEqIdx ] = std::max( maxCoeff[ contiFoamEqIdx ], std::abs( R2 ) / pvValue );
                }
                if constexpr (has_brine_) {
                    B_avg[ contiBrineEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiBrineEqIdx];
                    R_sum[ contiBrineEqIdx ] += R2;
                    maxCoeff[ contiBrineEqIdx ] = std::max( maxCoeff[ contiBrineEqIdx ], std::abs( R2 ) / pvValue );
                }

                if constexpr (has_polymermw_) {
                    static_assert(has_polymer_);

                    B_avg[contiPolymerMWEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    // the residual of the polymer molecular equation is scaled down by a 100, since molecular weight
                    // can be much bigger than 1, and this equation shares the same tolerance with other mass balance equations
                    // TODO: there should be a more general way to determine the scaling-down coefficient
                    const auto R2 = ebosResid[cell_idx][contiPolymerMWEqIdx] / 100.;
                    R_sum[contiPolymerMWEqIdx] += R2;
                    maxCoeff[contiPolymerMWEqIdx] = std::max( maxCoeff[contiPolymerMWEqIdx], std::abs( R2 ) / pvValue );
                }

                if constexpr (has_energy_) {
                    B_avg[ contiEnergyEqIdx ] += 1.0 / (4.182e1); // converting J -> RM3 (entalpy / (cp * deltaK * rho) assuming change of 1e-5K of water
                    const auto R2 = ebosResid[cell_idx][contiEnergyEqIdx];
                    R_sum[ contiEnergyEqIdx ] += R2;
                    maxCoeff[ contiEnergyEqIdx ] = std::max( maxCoeff[ contiEnergyEqIdx ], std::abs( R2 ) / pvValue );
                }

                if constexpr (has_micp_) {
                    B_avg[ contiMicrobialEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R1 = ebosResid[cell_idx][contiMicrobialEqIdx];
                    R_sum[ contiMicrobialEqIdx ] += R1;
                    maxCoeff[ contiMicrobialEqIdx ] = std::max( maxCoeff[ contiMicrobialEqIdx ], std::abs( R1 ) / pvValue );
                    B_avg[ contiOxygenEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiOxygenEqIdx];
                    R_sum[ contiOxygenEqIdx ] += R2;
                    maxCoeff[ contiOxygenEqIdx ] = std::max( maxCoeff[ contiOxygenEqIdx ], std::abs( R2 ) / pvValue );
                    B_avg[ contiUreaEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R3 = ebosResid[cell_idx][contiUreaEqIdx];
                    R_sum[ contiUreaEqIdx ] += R3;
                    maxCoeff[ contiUreaEqIdx ] = std::max( maxCoeff[ contiUreaEqIdx ], std::abs( R3 ) / pvValue );
                    B_avg[ contiBiofilmEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R4 = ebosResid[cell_idx][contiBiofilmEqIdx];
                    R_sum[ contiBiofilmEqIdx ] += R4;
                    maxCoeff[ contiBiofilmEqIdx ] = std::max( maxCoeff[ contiBiofilmEqIdx ], std::abs( R4 ) / pvValue );
                    B_avg[ contiCalciteEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R5 = ebosResid[cell_idx][contiCalciteEqIdx];
                    R_sum[ contiCalciteEqIdx ] += R5;
                    maxCoeff[ contiCalciteEqIdx ] = std::max( maxCoeff[ contiCalciteEqIdx ], std::abs( R5 ) / pvValue );
            }
          }

            OPM_END_PARALLEL_TRY_CATCH("BlackoilModelEbos::localConvergenceData() failed: ", grid_.comm());

            // compute local average in terms of global number of elements
            const int bSize = B_avg.size();
            for ( int i = 0; i<bSize; ++i )
            {
                B_avg[ i ] /= Scalar(domain.cells.size());
            }

            // return {pvSumLocal, numAquiferPvSumLocal};
            return pvSumLocal;
        }


        double computeCnvErrorPvLocal(const Domain& domain, const std::vector<Scalar>& B_avg, double dt)
        {
            double errorPV{};
            const auto& ebosModel = ebosSimulator_.model();
            const auto& ebosProblem = ebosSimulator_.problem();
            const auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            for (const int cell_idx : domain.cells)
            {
                const double pvValue = ebosProblem.referencePorosity(cell_idx, /*timeIdx=*/0) * ebosModel.dofTotalVolume( cell_idx );
                const auto& cellResidual = ebosResid[cell_idx];
                bool cnvViolated = false;

                for (unsigned eqIdx = 0; eqIdx < cellResidual.size(); ++eqIdx)
                {
                    using std::fabs;
                    Scalar CNV = cellResidual[eqIdx] * dt * B_avg[eqIdx] / pvValue;
                    cnvViolated = cnvViolated || (fabs(CNV) > param_.tolerance_cnv_);
                }

                if (cnvViolated)
                {
                    errorPV += pvValue;
                }
            }
            return errorPV;
        }


        /// \brief Compute the total pore volume of cells violating CNV that are not part
        ///        of a numerical aquifer.
        double computeCnvErrorPv(const std::vector<Scalar>& B_avg, double dt)
        {
            OPM_TIMEBLOCK(computeCnvErrorPv);
            double errorPV{};
            const auto& ebosModel = ebosSimulator_.model();
            const auto& ebosProblem = ebosSimulator_.problem();
            const auto& ebosResid = ebosSimulator_.model().linearizer().residual();
            const auto& gridView = ebosSimulator().gridView();
            ElementContext elemCtx(ebosSimulator_);

            OPM_BEGIN_PARALLEL_TRY_CATCH();

            for (const auto& elem : elements(gridView, Dune::Partitions::interiorBorder))
            {
                // Skip cells of numerical Aquifer
                if (isNumericalAquiferCell(gridView.grid(), elem))
                {
                    continue;
                }
                elemCtx.updatePrimaryStencil(elem);
                // elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const double pvValue = ebosProblem.referencePorosity(cell_idx, /*timeIdx=*/0) * ebosModel.dofTotalVolume( cell_idx );
                const auto& cellResidual = ebosResid[cell_idx];
                bool cnvViolated = false;

                for (unsigned eqIdx = 0; eqIdx < cellResidual.size(); ++eqIdx)
                {
                    using std::abs;
                    Scalar CNV = cellResidual[eqIdx] * dt * B_avg[eqIdx] / pvValue;
                    cnvViolated = cnvViolated || (abs(CNV) > param_.tolerance_cnv_);
                }

                if (cnvViolated)
                {
                    errorPV += pvValue;
                }
            }

            OPM_END_PARALLEL_TRY_CATCH("BlackoilModelEbos::ComputeCnvError() failed: ", grid_.comm());

            return grid_.comm().sum(errorPV);
        }

        ConvergenceReport getLocalReservoirConvergence(const double reportTime,
                                                       const double dt,
                                                       const int iteration,
                                                       const Domain& domain,
                                                       std::vector<Scalar>& B_avg,
                                                       std::vector<Scalar>& residual_norms)
        {
            typedef std::vector< Scalar > Vector;

            const int numComp = numEq;
            Vector R_sum(numComp, 0.0 );
            Vector maxCoeff(numComp, std::numeric_limits< Scalar >::lowest() );
            std::vector<int> maxCoeffCell(numComp, -1);
            const double pvSum = localDomainConvergenceData(domain, R_sum, maxCoeff, B_avg, maxCoeffCell);

            // compute global sum and max of quantities
            // const double pvSum = convergenceReduction(grid_.comm(), pvSumLocal,
            //                                           R_sum, maxCoeff, B_avg);

            auto cnvErrorPvFraction = computeCnvErrorPvLocal(domain, B_avg, dt);
            cnvErrorPvFraction /= pvSum;

            const double tol_mb  = param_.local_tolerance_scaling_mb_ * param_.tolerance_mb_;
            // Default value of relaxed_max_pv_fraction_ is 0.03 and min_strict_cnv_iter_ is 0.
            // For each iteration, we need to determine whether to use the relaxed CNV tolerance.
            // To disable the usage of relaxed CNV tolerance, you can set the relaxed_max_pv_fraction_ to be 0.
            const bool use_relaxed = cnvErrorPvFraction < param_.relaxed_max_pv_fraction_ && iteration >= param_.min_strict_cnv_iter_;
            // Tighter bound for local convergence should increase the
            // likelyhood of: local convergence => global convergence
            const double tol_cnv = param_.local_tolerance_scaling_cnv_
                * (use_relaxed ? param_.tolerance_cnv_relaxed_ :  param_.tolerance_cnv_);

            // Finish computation
            std::vector<Scalar> CNV(numComp);
            std::vector<Scalar> mass_balance_residual(numComp);
            for ( int compIdx = 0; compIdx < numComp; ++compIdx )
            {
                CNV[compIdx]                    = B_avg[compIdx] * dt * maxCoeff[compIdx];
                mass_balance_residual[compIdx]  = std::abs(B_avg[compIdx]*R_sum[compIdx]) * dt / pvSum;
                residual_norms.push_back(CNV[compIdx]);
            }

            // Create convergence report.
            ConvergenceReport report;
            using CR = ConvergenceReport;
            for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                double res[2] = { mass_balance_residual[compIdx], CNV[compIdx] };
                CR::ReservoirFailure::Type types[2] = { CR::ReservoirFailure::Type::MassBalance,
                                                        CR::ReservoirFailure::Type::Cnv };
                double tol[2] = { tol_mb, tol_cnv };
                int cell[2] = { -1, maxCoeffCell[compIdx] }; // No cell associated with MB failures.
                for (int ii : {0, 1}) {
                    if (std::isnan(res[ii])) {
                        report.setReservoirFailed({types[ii], CR::Severity::NotANumber, compIdx});//, cell[ii], res[ii]});
                        if ( terminal_output_ ) {
                            OpmLog::debug("NaN residual for " + compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] > maxResidualAllowed()) {
                        report.setReservoirFailed({types[ii], CR::Severity::TooLarge, compIdx});//, cell[ii], res[ii]});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Too large residual for " + compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] < 0.0) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});//, cell[ii], res[ii]});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Negative residual for " + compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] > tol[ii]) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});//, cell[ii], res[ii]});
                    }
                }
            }

            // Output of residuals.
            if ( terminal_output_ )
            {
                // Only rank 0 does print to std::cout
                if (iteration == 0) {
                    std::string msg = fmt::format("Domain {}, size {}, containing cell {}\n| Iter",
                                                  domain.index, domain.cells.size(), domain.cells[0]);
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    MB(";
                        msg += compNames_.name(compIdx)[0];
                        msg += ")  ";
                    }
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    CNV(";
                        msg += compNames_.name(compIdx)[0];
                        msg += ") ";
                    }
                    OpmLog::debug(msg);
                }
                std::ostringstream ss;
                ss << "| ";
                const std::streamsize oprec = ss.precision(3);
                const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                ss << std::setw(4) << iteration;
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << mass_balance_residual[compIdx];
                }
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << CNV[compIdx];
                }
                ss.precision(oprec);
                ss.flags(oflags);
                OpmLog::debug(ss.str());
            }

            return report;
        }



        ConvergenceReport getReservoirConvergence(const double reportTime,
                                                  const double dt,
                                                  const int iteration,
                                                  std::vector<Scalar>& B_avg,
                                                  std::vector<Scalar>& residual_norms)
        {
            OPM_TIMEBLOCK(getReservoirConvergence);
            typedef std::vector< Scalar > Vector;

            const int numComp = numEq;
            Vector R_sum(numComp, 0.0 );
            Vector maxCoeff(numComp, std::numeric_limits< Scalar >::lowest() );
            std::vector<int> maxCoeffCell(numComp, -1);
            const auto [ pvSumLocal, numAquiferPvSumLocal] = localConvergenceData(R_sum, maxCoeff, B_avg, maxCoeffCell);

            // compute global sum and max of quantities
            const auto [ pvSum, numAquiferPvSum ] =
                convergenceReduction(grid_.comm(), pvSumLocal,
                                     numAquiferPvSumLocal,
                                     R_sum, maxCoeff, B_avg);

            auto cnvErrorPvFraction = computeCnvErrorPv(B_avg, dt);
            cnvErrorPvFraction /= (pvSum - numAquiferPvSum);

            const double tol_mb  = param_.tolerance_mb_;
            // Default value of relaxed_max_pv_fraction_ is 0.03 and min_strict_cnv_iter_ is 0.
            // For each iteration, we need to determine whether to use the relaxed CNV tolerance.
            // To disable the usage of relaxed CNV tolerance, you can set the relaxed_max_pv_fraction_ to be 0.
            const bool use_relaxed = cnvErrorPvFraction < param_.relaxed_max_pv_fraction_ && iteration >= param_.min_strict_cnv_iter_;
            const double tol_cnv = use_relaxed ? param_.tolerance_cnv_relaxed_ :  param_.tolerance_cnv_;

            // Finish computation
            std::vector<Scalar> CNV(numComp);
            std::vector<Scalar> mass_balance_residual(numComp);
            for ( int compIdx = 0; compIdx < numComp; ++compIdx )
            {
                CNV[compIdx]                    = B_avg[compIdx] * dt * maxCoeff[compIdx];
                mass_balance_residual[compIdx]  = std::abs(B_avg[compIdx]*R_sum[compIdx]) * dt / pvSum;
                residual_norms.push_back(CNV[compIdx]);
            }

            // Create convergence report.
            ConvergenceReport report{reportTime};
            using CR = ConvergenceReport;
            for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                double res[2] = { mass_balance_residual[compIdx], CNV[compIdx] };
                CR::ReservoirFailure::Type types[2] = { CR::ReservoirFailure::Type::MassBalance,
                                                        CR::ReservoirFailure::Type::Cnv };
                double tol[2] = { tol_mb, tol_cnv };
                int cell[2] = { -1, maxCoeffCell[compIdx] }; // No cell associated with MB failures.
                for (int ii : {0, 1}) {
                    if (std::isnan(res[ii])) {
                        report.setReservoirFailed({types[ii], CR::Severity::NotANumber, compIdx});//, cell[ii], res[ii]});
                        if ( terminal_output_ ) {
                            OpmLog::debug("NaN residual for " + this->compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] > maxResidualAllowed()) {
                        report.setReservoirFailed({types[ii], CR::Severity::TooLarge, compIdx});//, cell[ii], res[ii]});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Too large residual for " + this->compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] < 0.0) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});//, cell[ii], res[ii]});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Negative residual for " + this->compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] > tol[ii]) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});//, cell[ii], res[ii]});
                    }
                    report.setReservoirConvergenceMetric(types[ii], compIdx, res[ii]);
                }
            }

            // Output of residuals.
            if ( terminal_output_ )
            {
                // Only rank 0 does print to std::cout
                if (iteration == 0) {
                    std::string msg = "Iter";
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    MB(";
                        msg += this->compNames_.name(compIdx)[0];
                        msg += ")  ";
                    }
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    CNV(";
                        msg += this->compNames_.name(compIdx)[0];
                        msg += ") ";
                    }
                    OpmLog::debug(msg);
                }
                std::ostringstream ss;
                const std::streamsize oprec = ss.precision(3);
                const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                ss << std::setw(4) << iteration;
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << mass_balance_residual[compIdx];
                }
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << CNV[compIdx];
                }
                ss.precision(oprec);
                ss.flags(oflags);
                OpmLog::debug(ss.str());
            }

            return report;
        }

        ConvergenceReport getLocalConvergence(const Domain& domain,
                                              const SimulatorTimerInterface& timer,
                                              const int iteration,
                                              std::vector<double>& residual_norms)
        {
            std::vector<Scalar> B_avg(numEq, 0.0);
            auto report = getLocalReservoirConvergence(timer.simulationTimeElapsed(),
                                                       timer.currentStepLength(),
                                                       iteration,
                                                       domain,
                                                       B_avg,
                                                       residual_norms);
            report += wellModel().getLocalWellConvergence(domain, B_avg, false);
            return report;
        }

        /// Compute convergence based on total mass balance (tol_mb) and maximum
        /// residual mass balance (tol_cnv).
        /// \param[in]   timer       simulation timer
        /// \param[in]   iteration   current iteration number
        /// \param[out]  residual_norms   CNV residuals by phase
        ConvergenceReport getConvergence(const SimulatorTimerInterface& timer,
                                         const int iteration,
                                         std::vector<double>& residual_norms)
        {
            OPM_TIMEBLOCK(getConvergence);
            // Get convergence reports for reservoir and wells.
            std::vector<Scalar> B_avg(numEq, 0.0);
            auto report = getReservoirConvergence(timer.simulationTimeElapsed(),
                                                  timer.currentStepLength(),
                                                  iteration, B_avg, residual_norms);
            {
                OPM_TIMEBLOCK(getWellConvergence);
                report += wellModel().getWellConvergence(B_avg, /*checkWellGroupControls*/report.converged());
            }
            return report;
        }


        /// The number of active fluid phases in the model.
        int numPhases() const
        {
            return phaseUsage_.num_phases;
        }

        /// Wrapper required due to not following generic API
        template<class T>
        std::vector<std::vector<double> >
        computeFluidInPlace(const T&, const std::vector<int>& fipnum) const
        {
            return computeFluidInPlace(fipnum);
        }

        /// Should not be called
        std::vector<std::vector<double> >
        computeFluidInPlace(const std::vector<int>& /*fipnum*/) const
        {
            OPM_TIMEBLOCK(computeFluidInPlace);
            //assert(true)
            //return an empty vector
            std::vector<std::vector<double> > regionValues(0, std::vector<double>(0,0.0));
            return regionValues;
        }

        const Simulator& ebosSimulator() const
        { return ebosSimulator_; }

        Simulator& ebosSimulator()
        { return ebosSimulator_; }

        /// return the statistics if the nonlinearIteration() method failed
        const SimulatorReportSingle& failureReport() const
        { return failureReport_; }

        /// return the statistics if the nonlinearIteration() method failed
        const SimulatorReportSingle& localAccumulatedReports() const
        { return local_reports_accumulated_; }

        const std::vector<StepReport>& stepReports() const
        {
            return convergence_reports_;
        }

        std::vector<StepReport> getStepReportsDestructively() const
        {
            return std::move(this->convergence_reports_);
        }

    protected:
        // ---------  Data members  ---------

        Simulator& ebosSimulator_;
        const Grid&            grid_;
        const PhaseUsage phaseUsage_;
        static constexpr bool has_solvent_ = getPropValue<TypeTag, Properties::EnableSolvent>();
        static constexpr bool has_extbo_ = getPropValue<TypeTag, Properties::EnableExtbo>();
        static constexpr bool has_polymer_ = getPropValue<TypeTag, Properties::EnablePolymer>();
        static constexpr bool has_polymermw_ = getPropValue<TypeTag, Properties::EnablePolymerMW>();
        static constexpr bool has_energy_ = getPropValue<TypeTag, Properties::EnableEnergy>();
        static constexpr bool has_foam_ = getPropValue<TypeTag, Properties::EnableFoam>();
        static constexpr bool has_brine_ = getPropValue<TypeTag, Properties::EnableBrine>();
        static constexpr bool has_micp_ = getPropValue<TypeTag, Properties::EnableMICP>();

        ModelParameters                 param_;
        SimulatorReportSingle failureReport_;
        SimulatorReportSingle local_reports_accumulated_;

        // Well Model
        BlackoilWellModel<TypeTag>& well_model_;

        /// \brief Whether we print something to std::cout
        bool terminal_output_;
        /// \brief The number of cells of the global grid.
        long int global_nc_;

        std::vector<std::vector<double>> residual_norms_history_;
        double current_relaxation_;
        BVector dx_old_;

        std::vector<StepReport> convergence_reports_;
        ComponentName compNames_{};
        std::vector<Domain> domains_;
        std::vector<std::unique_ptr<Mat>> domain_matrices_;
        std::vector<ISTLSolverType> domain_linsolvers_;

    public:
        /// return the StandardWells object
        BlackoilWellModel<TypeTag>&
        wellModel() { return well_model_; }

        const BlackoilWellModel<TypeTag>&
        wellModel() const { return well_model_; }

        void beginReportStep()
        {
            ebosSimulator_.problem().beginEpisode();
        }

        void endReportStep()
        {
            ebosSimulator_.problem().endEpisode();
        }

    private:
        template<class T>
        bool isNumericalAquiferCell(const Dune::CpGrid& grid, const T& elem)
        {
            const auto& aquiferCells = grid.sortedNumAquiferCells();
            if (aquiferCells.empty())
            {
                return false;
            }
            auto candidate = std::lower_bound(aquiferCells.begin(), aquiferCells.end(),
                                              elem.index());
            return candidate != aquiferCells.end() && *candidate == elem.index();
        }

        template<class G, class T>
        typename std::enable_if<!std::is_same<G,Dune::CpGrid>::value, bool>::type
        isNumericalAquiferCell(const G&, const T&)
        {
            return false;
        }

        double dpMaxRel() const { return param_.dp_max_rel_; }
        double dsMax() const { return param_.ds_max_; }
        double drMaxRel() const { return param_.dr_max_rel_; }
        double maxResidualAllowed() const { return param_.max_residual_allowed_; }
        double linear_solve_setup_time_;

    public:
        std::vector<bool> wasSwitched_;
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
