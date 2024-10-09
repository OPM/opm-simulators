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

#ifndef OPM_BLACKOILMODEL_NLDD_HEADER_INCLUDED
#define OPM_BLACKOILMODEL_NLDD_HEADER_INCLUDED

#include <dune/common/timer.hh>

#include <opm/grid/common/SubGridPart.hpp>

#include <opm/simulators/aquifers/AquiferGridUtils.hpp>

#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/flow/partitionCells.hpp>
#include <opm/simulators/flow/priVarsPacking.hpp>
#include <opm/simulators/flow/NonlinearSolver.hpp>
#include <opm/simulators/flow/SubDomain.hpp>

#include <opm/simulators/linalg/extractMatrix.hpp>

#if COMPILE_BDA_BRIDGE
#include <opm/simulators/linalg/ISTLSolverBda.hpp>
#else
#include <opm/simulators/linalg/ISTLSolver.hpp>
#endif

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>

#include <opm/simulators/utils/ComponentName.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm {

template<class TypeTag> class BlackoilModel;

/// A NLDD implementation for three-phase black oil.
template <class TypeTag>
class BlackoilModelNldd {
public:
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelParameters = BlackoilModelParameters<Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using BVector = typename BlackoilModel<TypeTag>::BVector;
    using Domain = SubDomain<Grid>;
    using ISTLSolverType = ISTLSolver<TypeTag>;
    using Mat = typename BlackoilModel<TypeTag>::Mat;

    static constexpr int numEq = Indices::numEq;

    //! \brief The constructor sets up the subdomains.
    //! \param model BlackOil model to solve for
    //! \param param param Model parameters
    //! \param compNames Names of the solution components
    BlackoilModelNldd(BlackoilModel<TypeTag>& model)
        : model_(model), rank_(model_.simulator().vanguard().grid().comm().rank())
    {
        // Create partitions.
        const auto& [partition_vector, num_domains] = this->partitionCells();

        // Scan through partitioning to get correct size for each.
        std::vector<int> sizes(num_domains, 0);
        for (const auto& p : partition_vector) {
            ++sizes[p];
        }

        // Set up correctly sized vectors of entity seeds and of indices for each partition.
        using EntitySeed = typename Grid::template Codim<0>::EntitySeed;
        std::vector<std::vector<EntitySeed>> seeds(num_domains);
        std::vector<std::vector<int>> partitions(num_domains);
        for (int domain = 0; domain < num_domains; ++domain) {
            seeds[domain].resize(sizes[domain]);
            partitions[domain].resize(sizes[domain]);
        }

        // Iterate through grid once, setting the seeds of all partitions.
        // Note: owned cells only!
        const auto& grid = model_.simulator().vanguard().grid();

        std::vector<int> count(num_domains, 0);
        const auto& gridView = grid.leafGridView();
        const auto beg = gridView.template begin<0, Dune::Interior_Partition>();
        const auto end = gridView.template end<0, Dune::Interior_Partition>();
        int cell = 0;
        for (auto it = beg; it != end; ++it, ++cell) {
            const int p = partition_vector[cell];
            seeds[p][count[p]] = it->seed();
            partitions[p][count[p]] = cell;
            ++count[p];
        }
        assert(count == sizes);

        // Create the domains.
        for (int index = 0; index < num_domains; ++index) {
            std::vector<bool> interior(partition_vector.size(), false);
            for (int ix : partitions[index]) {
                interior[ix] = true;
            }

            Dune::SubGridPart<Grid> view{grid, std::move(seeds[index])};

            this->domains_.emplace_back(index,
                                        std::move(partitions[index]),
                                        std::move(interior),
                                        std::move(view));
        }

        // Set up container for the local system matrices.
        domain_matrices_.resize(num_domains);

        // Set up container for the local linear solvers.
        for (int index = 0; index < num_domains; ++index) {
            // TODO: The ISTLSolver constructor will make
            // parallel structures appropriate for the full grid
            // only. This must be addressed before going parallel.
            const auto& eclState = model_.simulator().vanguard().eclState();
            FlowLinearSolverParameters loc_param;
            loc_param.is_nldd_local_solver_ = true;
            loc_param.init(eclState.getSimulationConfig().useCPR());
            // Override solver type with umfpack if small domain.
            if (domains_[index].cells.size() < 200) {
                loc_param.linsolver_ = "umfpack";
            }
            loc_param.linear_solver_print_json_definition_ = false;
            const bool force_serial = true;
            domain_linsolvers_.emplace_back(model_.simulator(), loc_param, force_serial);
        }

        assert(int(domains_.size()) == num_domains);
    }

    //! \brief Called before starting a time step.
    void prepareStep()
    {
        // Setup domain->well mapping.
        model_.wellModel().setupDomains(domains_);
    }

    //! \brief Do one non-linear NLDD iteration.
    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationNldd(const int iteration,
                                                 const SimulatorTimerInterface& timer,
                                                 NonlinearSolverType& nonlinear_solver)
    {
        // -----------   Set up reports and timer   -----------
        SimulatorReportSingle report;
        Dune::Timer perfTimer;

        model_.initialLinearization(report, iteration, nonlinear_solver.minIter(), nonlinear_solver.maxIter(), timer);

        if (report.converged) {
            return report;
        }

        // -----------   If not converged, do an NLDD iteration   -----------

        auto& solution = model_.simulator().model().solution(0);
        auto initial_solution = solution;
        auto locally_solved = initial_solution;

        // -----------   Decide on an ordering for the domains   -----------
        const auto domain_order = this->getSubdomainOrder();

        // -----------   Solve each domain separately   -----------
        DeferredLogger logger;
        std::vector<SimulatorReportSingle> domain_reports(domains_.size());
        for (const int domain_index : domain_order) {
            const auto& domain = domains_[domain_index];
            SimulatorReportSingle local_report;
            try {
                switch (model_.param().local_solve_approach_) {
                case DomainSolveApproach::Jacobi:
                    solveDomainJacobi(solution, locally_solved, local_report, logger,
                                      iteration, timer, domain);
                    break;
                default:
                case DomainSolveApproach::GaussSeidel:
                    solveDomainGaussSeidel(solution, locally_solved, local_report, logger,
                                           iteration, timer, domain);
                    break;
                }
            }
            catch (...) {
                // Something went wrong during local solves.
                local_report.converged = false;
            }
            // This should have updated the global matrix to be
            // dR_i/du_j evaluated at new local solutions for
            // i == j, at old solution for i != j.
            if (!local_report.converged) {
                // TODO: more proper treatment, including in parallel.
                logger.debug(fmt::format("Convergence failure in domain {} on rank {}." , domain.index, rank_));
            }
            domain_reports[domain.index] = local_report;
        }

        // Communicate and log all messages.
        auto global_logger = gatherDeferredLogger(logger, model_.simulator().vanguard().grid().comm());
        global_logger.logMessages();

        // Accumulate local solve data.
        // Putting the counts in a single array to avoid multiple
        // comm.sum() calls. Keeping the named vars for readability.
        std::array<int, 4> counts{ 0, 0, 0, static_cast<int>(domain_reports.size()) };
        int& num_converged = counts[0];
        int& num_converged_already = counts[1];
        int& num_local_newtons = counts[2];
        int& num_domains = counts[3];
        {
            SimulatorReportSingle rep;
            for (const auto& dr : domain_reports) {
                if (dr.converged) {
                    ++num_converged;
                    if (dr.total_newton_iterations == 0) {
                        ++num_converged_already;
                    }
                }
                rep += dr;
            }
            num_local_newtons = rep.total_newton_iterations;
            local_reports_accumulated_ += rep;
        }

        if (model_.param().local_solve_approach_ == DomainSolveApproach::Jacobi) {
            solution = locally_solved;
            model_.simulator().model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }

#if HAVE_MPI
        // Communicate solutions:
        // With multiple processes, this process' overlap (i.e. not
        // owned) cells' solution values have been modified by local
        // solves in the owning processes, and remain unchanged
        // here. We must therefore receive the updated solution on the
        // overlap cells and update their intensive quantities before
        // we move on.
        const auto& comm = model_.simulator().vanguard().grid().comm();
        if (comm.size() > 1) {
            const auto* ccomm = model_.simulator().model().newtonMethod().linearSolver().comm();

            // Copy numerical values from primary vars.
            ccomm->copyOwnerToAll(solution, solution);

            // Copy flags from primary vars.
            const std::size_t num = solution.size();
            Dune::BlockVector<std::size_t> allmeanings(num);
            for (std::size_t ii = 0; ii < num; ++ii) {
                allmeanings[ii] = PVUtil::pack(solution[ii]);
            }
            ccomm->copyOwnerToAll(allmeanings, allmeanings);
            for (std::size_t ii = 0; ii < num; ++ii) {
                PVUtil::unPack(solution[ii], allmeanings[ii]);
            }

            // Update intensive quantities for our overlap values.
            model_.simulator().model().invalidateAndUpdateIntensiveQuantitiesOverlap(/*timeIdx=*/0);

            // Make total counts of domains converged.
            comm.sum(counts.data(), counts.size());
        }
#endif // HAVE_MPI

        const bool is_iorank = this->rank_ == 0;
        if (is_iorank) {
            OpmLog::debug(fmt::format("Local solves finished. Converged for {}/{} domains. {} domains did no work. {} total local Newton iterations.\n",
                                      num_converged, num_domains, num_converged_already, num_local_newtons));
        }

        // Finish with a Newton step.
        // Note that the "iteration + 100" is a simple way to avoid entering
        // "if (iteration == 0)" and similar blocks, and also makes it a little
        // easier to spot the iteration residuals in the DBG file. A more sophisticated
        // approach can be done later.
        auto rep = model_.nonlinearIterationNewton(iteration + 100, timer, nonlinear_solver);
        report += rep;
        if (rep.converged) {
            report.converged = true;
        }
        return report;
    }

    /// return the statistics if the nonlinearIteration() method failed
    const SimulatorReportSingle& localAccumulatedReports() const
    {
        return local_reports_accumulated_;
    }

    void writePartitions(const std::filesystem::path& odir) const
    {
        const auto& elementMapper = this->model_.simulator().model().elementMapper();
        const auto& cartMapper = this->model_.simulator().vanguard().cartesianIndexMapper();

        const auto& grid = this->model_.simulator().vanguard().grid();
        const auto& comm = grid.comm();
        const auto nDigit = 1 + static_cast<int>(std::floor(std::log10(comm.size())));

        std::ofstream pfile { odir / fmt::format("{1:0>{0}}", nDigit, comm.rank()) };

        const auto p = this->reconstitutePartitionVector();
        auto i = 0;
        for (const auto& cell : elements(grid.leafGridView(), Dune::Partitions::interior)) {
            pfile << comm.rank() << ' '
                  << cartMapper.cartesianIndex(elementMapper.index(cell)) << ' '
                  << p[i++] << '\n';
        }
    }

private:

    //! \brief Solve the equation system for a single domain.
    std::pair<SimulatorReportSingle, ConvergenceReport>
    solveDomain(const Domain& domain,
                const SimulatorTimerInterface& timer,
                DeferredLogger& logger,
                [[maybe_unused]] const int global_iteration,
                const bool initial_assembly_required)
    {
        auto& modelSimulator = model_.simulator();

        SimulatorReportSingle report;
        Dune::Timer solveTimer;
        solveTimer.start();
        Dune::Timer detailTimer;

        modelSimulator.model().newtonMethod().setIterationIndex(0);

        // When called, if assembly has already been performed
        // with the initial values, we only need to check
        // for local convergence. Otherwise, we must do a local
        // assembly.
        int iter = 0;
        if (initial_assembly_required) {
            detailTimer.start();
            modelSimulator.model().newtonMethod().setIterationIndex(iter);
            // TODO: we should have a beginIterationLocal function()
            // only handling the well model for now
            modelSimulator.problem().wellModel().assembleDomain(modelSimulator.model().newtonMethod().numIterations(),
                                                                modelSimulator.timeStepSize(),
                                                               domain);
            // Assemble reservoir locally.
            report += this->assembleReservoirDomain(domain);
            report.assemble_time += detailTimer.stop();
        }
        detailTimer.reset();
        detailTimer.start();
        std::vector<Scalar> resnorms;
        auto convreport = this->getDomainConvergence(domain, timer, 0, logger, resnorms);
        if (convreport.converged()) {
            // TODO: set more info, timing etc.
            report.converged = true;
            return { report, convreport };
        }

        // We have already assembled for the first iteration,
        // but not done the Schur complement for the wells yet.
        detailTimer.reset();
        detailTimer.start();
        model_.wellModel().linearizeDomain(domain,
                                           modelSimulator.model().linearizer().jacobian(),
                                           modelSimulator.model().linearizer().residual());
        const double tt1 = detailTimer.stop();
        report.assemble_time += tt1;
        report.assemble_time_well += tt1;

        // Local Newton loop.
        const int max_iter = model_.param().max_local_solve_iterations_;
        const auto& grid = modelSimulator.vanguard().grid();
        double damping_factor = 1.0;
        std::vector<std::vector<Scalar>> convergence_history;
        convergence_history.reserve(20);
        convergence_history.push_back(resnorms);
        do {
            // Solve local linear system.
            // Note that x has full size, we expect it to be nonzero only for in-domain cells.
            const int nc = grid.size(0);
            BVector x(nc);
            detailTimer.reset();
            detailTimer.start();
            this->solveJacobianSystemDomain(domain, x);
            model_.wellModel().postSolveDomain(x, domain);
            if (damping_factor != 1.0) {
                x *= damping_factor;
            }
            report.linear_solve_time += detailTimer.stop();
            report.linear_solve_setup_time += model_.linearSolveSetupTime();
            report.total_linear_iterations = model_.linearIterationsLastSolve();

            // Update local solution. // TODO: x is still full size, should we optimize it?
            detailTimer.reset();
            detailTimer.start();
            this->updateDomainSolution(domain, x);
            report.update_time += detailTimer.stop();

            // Assemble well and reservoir.
            detailTimer.reset();
            detailTimer.start();
            ++iter;
            modelSimulator.model().newtonMethod().setIterationIndex(iter);
            // TODO: we should have a beginIterationLocal function()
            // only handling the well model for now
            // Assemble reservoir locally.
            modelSimulator.problem().wellModel().assembleDomain(modelSimulator.model().newtonMethod().numIterations(),
                                                                modelSimulator.timeStepSize(),
                                                                domain);
            report += this->assembleReservoirDomain(domain);
            report.assemble_time += detailTimer.stop();

            // Check for local convergence.
            detailTimer.reset();
            detailTimer.start();
            resnorms.clear();
            convreport = this->getDomainConvergence(domain, timer, iter, logger, resnorms);
            convergence_history.push_back(resnorms);

            // apply the Schur complement of the well model to the
            // reservoir linearized equations
            detailTimer.reset();
            detailTimer.start();
            model_.wellModel().linearizeDomain(domain,
                                               modelSimulator.model().linearizer().jacobian(),
                                               modelSimulator.model().linearizer().residual());
            const double tt2 = detailTimer.stop();
            report.assemble_time += tt2;
            report.assemble_time_well += tt2;

            // Check if we should dampen. Only do so if wells are converged.
            if (!convreport.converged() && !convreport.wellFailed()) {
                bool oscillate = false;
                bool stagnate = false;
                const int numPhases = convergence_history.front().size();
                detail::detectOscillations(convergence_history, iter, numPhases,
                                           Scalar{0.2}, 1, oscillate, stagnate);
                if (oscillate) {
                    damping_factor *= 0.85;
                    logger.debug(fmt::format("| Damping factor is now {}", damping_factor));
                }
            }
        } while (!convreport.converged() && iter <= max_iter);

        modelSimulator.problem().endIteration();

        report.converged = convreport.converged();
        report.total_newton_iterations = iter;
        report.total_linearizations = iter;
        report.total_time = solveTimer.stop();
        // TODO: set more info, timing etc.
        return { report, convreport };
    }

    /// Assemble the residual and Jacobian of the nonlinear system.
    SimulatorReportSingle assembleReservoirDomain(const Domain& domain)
    {
        // -------- Mass balance equations --------
        model_.simulator().model().linearizer().linearizeDomain(domain);
        return model_.wellModel().lastReport();
    }

    //! \brief Solve the linearized system for a domain.
    void solveJacobianSystemDomain(const Domain& domain, BVector& global_x)
    {
        const auto& modelSimulator = model_.simulator();

        Dune::Timer perfTimer;
        perfTimer.start();

        const Mat& main_matrix = modelSimulator.model().linearizer().jacobian().istlMatrix();
        if (domain_matrices_[domain.index]) {
            Details::copySubMatrix(main_matrix, domain.cells, *domain_matrices_[domain.index]);
        } else {
            domain_matrices_[domain.index] = std::make_unique<Mat>(Details::extractMatrix(main_matrix, domain.cells));
        }
        auto& jac = *domain_matrices_[domain.index];
        auto res = Details::extractVector(modelSimulator.model().linearizer().residual(),
                                          domain.cells);
        auto x = res;

        // set initial guess
        global_x = 0.0;
        x = 0.0;

        auto& linsolver = domain_linsolvers_[domain.index];

        linsolver.prepare(jac, res);
        model_.linearSolveSetupTime() = perfTimer.stop();
        linsolver.setResidual(res);
        linsolver.solve(x);

        Details::setGlobal(x, domain.cells, global_x);
    }

    /// Apply an update to the primary variables.
    void updateDomainSolution(const Domain& domain, const BVector& dx)
    {
        auto& simulator = model_.simulator();
        auto& newtonMethod = simulator.model().newtonMethod();
        SolutionVector& solution = simulator.model().solution(/*timeIdx=*/0);

        newtonMethod.update_(/*nextSolution=*/solution,
                             /*curSolution=*/solution,
                             /*update=*/dx,
                             /*resid=*/dx,
                             domain.cells); // the update routines of the black
                                            // oil model do not care about the
                                            // residual

        // if the solution is updated, the intensive quantities need to be recalculated
        simulator.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain);
    }

    //! \brief Get reservoir quantities on this process needed for convergence calculations.
    std::pair<Scalar, Scalar> localDomainConvergenceData(const Domain& domain,
                                                         std::vector<Scalar>& R_sum,
                                                         std::vector<Scalar>& maxCoeff,
                                                         std::vector<Scalar>& B_avg,
                                                         std::vector<int>& maxCoeffCell)
    {
        const auto& modelSimulator = model_.simulator();

        Scalar pvSumLocal = 0.0;
        Scalar numAquiferPvSumLocal = 0.0;
        const auto& model = modelSimulator.model();
        const auto& problem = modelSimulator.problem();

        const auto& modelResid = modelSimulator.model().linearizer().residual();

        ElementContext elemCtx(modelSimulator);
        const auto& gridView = domain.view;
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        IsNumericalAquiferCell isNumericalAquiferCell(gridView.grid());

        for (auto elemIt = gridView.template begin</*codim=*/0>();
             elemIt != elemEndIt;
             ++elemIt)
        {
            if (elemIt->partitionType() != Dune::InteriorEntity) {
                continue;
            }
            const auto& elem = *elemIt;
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            const auto pvValue = problem.referencePorosity(cell_idx, /*timeIdx=*/0) *
                                 model.dofTotalVolume(cell_idx);
            pvSumLocal += pvValue;

            if (isNumericalAquiferCell(elem))
            {
                numAquiferPvSumLocal += pvValue;
            }

            model_.getMaxCoeff(cell_idx, intQuants, fs, modelResid, pvValue,
                               B_avg, R_sum, maxCoeff, maxCoeffCell);
        }

        // compute local average in terms of global number of elements
        const int bSize = B_avg.size();
        for ( int i = 0; i<bSize; ++i )
        {
            B_avg[ i ] /= Scalar(domain.cells.size());
        }

        return {pvSumLocal, numAquiferPvSumLocal};
    }

    ConvergenceReport getDomainReservoirConvergence(const double reportTime,
                                                    const double dt,
                                                    const int iteration,
                                                    const Domain& domain,
                                                    DeferredLogger& logger,
                                                    std::vector<Scalar>& B_avg,
                                                    std::vector<Scalar>& residual_norms)
    {
        using Vector = std::vector<Scalar>;

        const int numComp = numEq;
        Vector R_sum(numComp, 0.0 );
        Vector maxCoeff(numComp, std::numeric_limits<Scalar>::lowest() );
        std::vector<int> maxCoeffCell(numComp, -1);
        const auto [ pvSum, numAquiferPvSum]
            = this->localDomainConvergenceData(domain, R_sum, maxCoeff, B_avg, maxCoeffCell);

        auto cnvErrorPvFraction = computeCnvErrorPvLocal(domain, B_avg, dt);
        cnvErrorPvFraction /= (pvSum - numAquiferPvSum);

        // Default value of relaxed_max_pv_fraction_ is 0.03 and min_strict_cnv_iter_ is 0.
        // For each iteration, we need to determine whether to use the relaxed CNV tolerance.
        // To disable the usage of relaxed CNV tolerance, you can set the relaxed_max_pv_fraction_ to be 0.
        const bool use_relaxed_cnv = cnvErrorPvFraction < model_.param().relaxed_max_pv_fraction_ &&
                                 iteration >= model_.param().min_strict_cnv_iter_;
        // Tighter bound for local convergence should increase the
        // likelyhood of: local convergence => global convergence
        const Scalar tol_cnv = model_.param().local_tolerance_scaling_cnv_
            * (use_relaxed_cnv ? model_.param().tolerance_cnv_relaxed_
                           : model_.param().tolerance_cnv_);

        const bool use_relaxed_mb = iteration >= model_.param().min_strict_mb_iter_;
        const Scalar tol_mb  = model_.param().local_tolerance_scaling_mb_
            * (use_relaxed_mb ? model_.param().tolerance_mb_relaxed_ : model_.param().tolerance_mb_);

        // Finish computation
        std::vector<Scalar> CNV(numComp);
        std::vector<Scalar> mass_balance_residual(numComp);
        for (int compIdx = 0; compIdx < numComp; ++compIdx )
        {
            CNV[compIdx]                    = B_avg[compIdx] * dt * maxCoeff[compIdx];
            mass_balance_residual[compIdx]  = std::abs(B_avg[compIdx]*R_sum[compIdx]) * dt / pvSum;
            residual_norms.push_back(CNV[compIdx]);
        }

        // Create convergence report.
        ConvergenceReport report{reportTime};
        using CR = ConvergenceReport;
        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            Scalar res[2] = { mass_balance_residual[compIdx], CNV[compIdx] };
            CR::ReservoirFailure::Type types[2] = { CR::ReservoirFailure::Type::MassBalance,
                                                    CR::ReservoirFailure::Type::Cnv };
            Scalar tol[2] = { tol_mb, tol_cnv };
            for (int ii : {0, 1}) {
                if (std::isnan(res[ii])) {
                    report.setReservoirFailed({types[ii], CR::Severity::NotANumber, compIdx});
                    logger.debug("NaN residual for " + model_.compNames().name(compIdx) + " equation.");
                } else if (res[ii] > model_.param().max_residual_allowed_) {
                    report.setReservoirFailed({types[ii], CR::Severity::TooLarge, compIdx});
                    logger.debug("Too large residual for " + model_.compNames().name(compIdx) + " equation.");
                } else if (res[ii] < 0.0) {
                    report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                    logger.debug("Negative residual for " + model_.compNames().name(compIdx) + " equation.");
                } else if (res[ii] > tol[ii]) {
                    report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                }

                report.setReservoirConvergenceMetric(types[ii], compIdx, res[ii], tol[ii]);
            }
        }

        // Output of residuals. If converged at initial state, log nothing.
        const bool converged_at_initial_state = (report.converged() && iteration == 0);
        if (!converged_at_initial_state) {
            if (iteration == 0) {
                // Log header.
                std::string msg = fmt::format("Domain {} on rank {}, size {}, containing cell {}\n| Iter",
                                              domain.index, this->rank_, domain.cells.size(), domain.cells[0]);
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    msg += "    MB(";
                    msg += model_.compNames().name(compIdx)[0];
                    msg += ")  ";
                }
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    msg += "    CNV(";
                    msg += model_.compNames().name(compIdx)[0];
                    msg += ") ";
                }
                logger.debug(msg);
            }
            // Log convergence data.
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
            logger.debug(ss.str());
        }

        return report;
    }

    ConvergenceReport getDomainConvergence(const Domain& domain,
                                           const SimulatorTimerInterface& timer,
                                           const int iteration,
                                           DeferredLogger& logger,
                                           std::vector<Scalar>& residual_norms)
    {
        std::vector<Scalar> B_avg(numEq, 0.0);
        auto report = this->getDomainReservoirConvergence(timer.simulationTimeElapsed(),
                                                          timer.currentStepLength(),
                                                          iteration,
                                                          domain,
                                                          logger,
                                                          B_avg,
                                                          residual_norms);
        report += model_.wellModel().getDomainWellConvergence(domain, B_avg, logger);
        return report;
    }

    //! \brief Returns subdomain ordered according to method and ordering measure.
    std::vector<int> getSubdomainOrder()
    {
        const auto& modelSimulator = model_.simulator();
        const auto& solution = modelSimulator.model().solution(0);

        std::vector<int> domain_order(domains_.size());
        std::iota(domain_order.begin(), domain_order.end(), 0);

        if (model_.param().local_solve_approach_ == DomainSolveApproach::Jacobi) {
            // Do nothing, 0..n-1 order is fine.
            return domain_order;
        } else if (model_.param().local_solve_approach_ == DomainSolveApproach::GaussSeidel) {
            // Calculate the measure used to order the domains.
            std::vector<Scalar> measure_per_domain(domains_.size());
            switch (model_.param().local_domain_ordering_) {
            case DomainOrderingMeasure::AveragePressure: {
                // Use average pressures to order domains.
                for (const auto& domain : domains_) {
                    Scalar press_sum = 0.0;
                    for (const int c : domain.cells) {
                        press_sum += solution[c][Indices::pressureSwitchIdx];
                    }
                    const Scalar avgpress = press_sum / domain.cells.size();
                    measure_per_domain[domain.index] = avgpress;
                }
                break;
            }
            case DomainOrderingMeasure::MaxPressure: {
                // Use max pressures to order domains.
                for (const auto& domain : domains_) {
                    Scalar maxpress = 0.0;
                    for (const int c : domain.cells) {
                        maxpress = std::max(maxpress, solution[c][Indices::pressureSwitchIdx]);
                    }
                    measure_per_domain[domain.index] = maxpress;
                }
                break;
            }
            case DomainOrderingMeasure::Residual: {
                // Use maximum residual to order domains.
                const auto& residual = modelSimulator.model().linearizer().residual();
                const int num_vars = residual[0].size();
                for (const auto& domain : domains_) {
                    Scalar maxres = 0.0;
                    for (const int c : domain.cells) {
                        for (int ii = 0; ii < num_vars; ++ii) {
                            maxres = std::max(maxres, std::fabs(residual[c][ii]));
                        }
                    }
                    measure_per_domain[domain.index] = maxres;
                }
                break;
            }
            } // end of switch (model_.param().local_domain_ordering_)

            // Sort by largest measure, keeping index order if equal.
            const auto& m = measure_per_domain;
            std::stable_sort(domain_order.begin(), domain_order.end(),
                             [&m](const int i1, const int i2){ return m[i1] > m[i2]; });
            return domain_order;
        } else {
            throw std::logic_error("Domain solve approach must be Jacobi or Gauss-Seidel");
        }
    }

    template<class GlobalEqVector>
    void solveDomainJacobi(GlobalEqVector& solution,
                           GlobalEqVector& locally_solved,
                           SimulatorReportSingle& local_report,
                           DeferredLogger& logger,
                           const int iteration,
                           const SimulatorTimerInterface& timer,
                           const Domain& domain)
    {
        auto initial_local_well_primary_vars = model_.wellModel().getPrimaryVarsDomain(domain);
        auto initial_local_solution = Details::extractVector(solution, domain.cells);
        auto res = solveDomain(domain, timer, logger, iteration, false);
        local_report = res.first;
        if (local_report.converged) {
            auto local_solution = Details::extractVector(solution, domain.cells);
            Details::setGlobal(local_solution, domain.cells, locally_solved);
            Details::setGlobal(initial_local_solution, domain.cells, solution);
            model_.simulator().model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain);
        } else {
            model_.wellModel().setPrimaryVarsDomain(domain, initial_local_well_primary_vars);
            Details::setGlobal(initial_local_solution, domain.cells, solution);
            model_.simulator().model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain);
        }
    }

    template<class GlobalEqVector>
    void solveDomainGaussSeidel(GlobalEqVector& solution,
                                GlobalEqVector& locally_solved,
                                SimulatorReportSingle& local_report,
                                DeferredLogger& logger,
                                const int iteration,
                                const SimulatorTimerInterface& timer,
                                const Domain& domain)
    {
        auto initial_local_well_primary_vars = model_.wellModel().getPrimaryVarsDomain(domain);
        auto initial_local_solution = Details::extractVector(solution, domain.cells);
        auto res = solveDomain(domain, timer, logger, iteration, true);
        local_report = res.first;
        if (!local_report.converged) {
            // We look at the detailed convergence report to evaluate
            // if we should accept the unconverged solution.
            const auto& convrep = res.second;
            // We do not accept a solution if the wells are unconverged.
            if (!convrep.wellFailed()) {
                // Calculare the sums of the mb and cnv failures.
                Scalar mb_sum = 0.0;
                Scalar cnv_sum = 0.0;
                for (const auto& rc : convrep.reservoirConvergence()) {
                    if (rc.type() == ConvergenceReport::ReservoirFailure::Type::MassBalance) {
                        mb_sum += rc.value();
                    } else if (rc.type() == ConvergenceReport::ReservoirFailure::Type::Cnv) {
                        cnv_sum += rc.value();
                    }
                }
                // If not too high, we overrule the convergence failure.
                const Scalar acceptable_local_mb_sum = 1e-3;
                const Scalar acceptable_local_cnv_sum = 1.0;
                if (mb_sum < acceptable_local_mb_sum && cnv_sum < acceptable_local_cnv_sum) {
                    local_report.converged = true;
                    logger.debug(fmt::format("Accepting solution in unconverged domain {} on rank {}.", domain.index, rank_));
                    logger.debug(fmt::format("Value of mb_sum: {}  cnv_sum: {}", mb_sum, cnv_sum));
                } else {
                    logger.debug("Unconverged local solution.");
                }
            } else {
                logger.debug("Unconverged local solution with well convergence failures:");
                for (const auto& wf : convrep.wellFailures()) {
                    logger.debug(to_string(wf));
                }
            }
        }
        if (local_report.converged) {
            auto local_solution = Details::extractVector(solution, domain.cells);
            Details::setGlobal(local_solution, domain.cells, locally_solved);
        } else {
            model_.wellModel().setPrimaryVarsDomain(domain, initial_local_well_primary_vars);
            Details::setGlobal(initial_local_solution, domain.cells, solution);
            model_.simulator().model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0, domain);
        }
    }

    Scalar computeCnvErrorPvLocal(const Domain& domain,
                                  const std::vector<Scalar>& B_avg, double dt) const
    {
        Scalar errorPV{};
        const auto& simulator = model_.simulator();
        const auto& model = simulator.model();
        const auto& problem = simulator.problem();
        const auto& residual = simulator.model().linearizer().residual();

        for (const int cell_idx : domain.cells) {
            const Scalar pvValue = problem.referencePorosity(cell_idx, /*timeIdx=*/0) *
                                   model.dofTotalVolume(cell_idx);
            const auto& cellResidual = residual[cell_idx];
            bool cnvViolated = false;

            for (unsigned eqIdx = 0; eqIdx < cellResidual.size(); ++eqIdx) {
                using std::fabs;
                Scalar CNV = cellResidual[eqIdx] * dt * B_avg[eqIdx] / pvValue;
                cnvViolated = cnvViolated || (fabs(CNV) > model_.param().tolerance_cnv_);
            }

            if (cnvViolated) {
                errorPV += pvValue;
            }
        }
        return errorPV;
    }

    decltype(auto) partitionCells() const
    {
        const auto& grid = this->model_.simulator().vanguard().grid();

        using GridView = std::remove_cv_t<std::remove_reference_t<
            decltype(grid.leafGridView())>>;

        using Element = std::remove_cv_t<std::remove_reference_t<
            typename GridView::template Codim<0>::Entity>>;

        const auto& param = this->model_.param();

        auto zoltan_ctrl = ZoltanPartitioningControl<Element>{};

        zoltan_ctrl.domain_imbalance = param.local_domain_partition_imbalance_;

        zoltan_ctrl.index =
            [elementMapper = &this->model_.simulator().model().elementMapper()]
            (const Element& element)
        {
            return elementMapper->index(element);
        };

        zoltan_ctrl.local_to_global =
            [cartMapper = &this->model_.simulator().vanguard().cartesianIndexMapper()]
            (const int elemIdx)
        {
            return cartMapper->cartesianIndex(elemIdx);
        };

        // Forming the list of wells is expensive, so do this only if needed.
        const auto need_wells = param.local_domain_partition_method_ == "zoltan";

        const auto wells = need_wells
            ? this->model_.simulator().vanguard().schedule().getWellsatEnd()
            : std::vector<Well>{};

        const auto& possibleFutureConnectionSet = need_wells
            ? this->model_.simulator().vanguard().schedule().getPossibleFutureConnections()
            : std::unordered_map<std::string, std::set<int>> {};

        // If defaulted parameter for number of domains, choose a reasonable default.
        constexpr int default_cells_per_domain = 1000;
        const int num_domains = (param.num_local_domains_ > 0)
            ? param.num_local_domains_
            : detail::countGlobalCells(grid) / default_cells_per_domain;

        return ::Opm::partitionCells(param.local_domain_partition_method_,
                                     num_domains, grid.leafGridView(), wells,
                                     possibleFutureConnectionSet, zoltan_ctrl);
    }

    std::vector<int> reconstitutePartitionVector() const
    {
        const auto& grid = this->model_.simulator().vanguard().grid();

        auto numD = std::vector<int>(grid.comm().size() + 1, 0);
        numD[grid.comm().rank() + 1] = static_cast<int>(this->domains_.size());
        grid.comm().sum(numD.data(), numD.size());
        std::partial_sum(numD.begin(), numD.end(), numD.begin());

        auto p = std::vector<int>(grid.size(0));
        auto maxCellIdx = std::numeric_limits<int>::min();

        auto d = numD[grid.comm().rank()];
        for (const auto& domain : this->domains_) {
            for (const auto& cell : domain.cells) {
                p[cell] = d;
                if (cell > maxCellIdx) {
                    maxCellIdx = cell;
                }
            }

            ++d;
        }

        p.erase(p.begin() + maxCellIdx + 1, p.end());
        return p;
    }

    BlackoilModel<TypeTag>& model_; //!< Reference to model
    std::vector<Domain> domains_; //!< Vector of subdomains
    std::vector<std::unique_ptr<Mat>> domain_matrices_; //!< Vector of matrix operator for each subdomain
    std::vector<ISTLSolverType> domain_linsolvers_; //!< Vector of linear solvers for each domain
    SimulatorReportSingle local_reports_accumulated_; //!< Accumulated convergence report for subdomain solvers
    int rank_ = 0; //!< MPI rank of this process
};

} // namespace Opm

#endif // OPM_BLACKOILMODEL_NLDD_HEADER_INCLUDED
