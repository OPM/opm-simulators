/*
  Copyright 2016 - 2019 SINTEF Digital, Mathematics & Cybernetics.
  Copyright 2016 - 2018 Equinor ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 Norce AS

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

#ifndef OPM_BLACKOILWELLMODEL_NLDD_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_NLDD_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_BLACKOILWELLMODEL_NLDD_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelNldd.hpp>
#endif

#include <algorithm>

namespace Opm {

template<typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
assemble(const double dt,
         const Domain& domain)
{
    OPM_TIMEBLOCK(assemble);
    // NLDD domain well assembly for local solves.
    //
    // Call chain:  BlackoilWellModelNldd::assemble()
    //                -> BlackoilWellModelNldd::assembleWellEq(dt, domain)
    //                   -> WellInterface::assembleWellEq()
    //                      -> WellInterface::prepareWellBeforeAssembling()
    //                      -> WellInterface::assembleWellEqWithoutIteration()
    //
    // Key assumptions and design decisions:
    //
    // 1. Global well initialization (BlackoilWellModel::prepareTimeStep, including the
    //    optional initial well solve) has already been called via
    //    BlackoilWellModel::assemble() during the initial global linearization.
    //    This is ensured by NewtonIterationContext::needsTimestepInit().
    //
    // 2. Inner well iterations (BlackoilWellModel::iterateWellEquations) are skipped
    //    during local solves via NewtonIterationContext::shouldRunInnerWellIterations()
    //    (inLocalSolve=true).
    //    The NLDD local Newton loop already iterates wells to convergence.
    //
    // 3. No MPI collective operations (do_mpi_gather=false).
    //
    // 4. Only individual well constraints are checked (see
    //    BlackoilWellModelNldd::updateWellControls).
    //    Group controls, network balancing, gas lift, and NUPCOL/VREP updates
    //    are NOT performed during domain solves. This means:
    //    - Wells can switch between individual controls (BHP, THP, rate limits)
    //      to respect physical safety constraints.
    //    - Wells cannot be switched to/from GRUP control.  Individual constraint
    //      checks (WellInterface::checkIndividualConstraints) never produce GRUP
    //      as a result; only WellInterface::checkGroupConstraints can do that,
    //      and it is not called here.
    //    - Group production targets and VREP/REIN injection targets use stale
    //      group state from the last global assembly.
    //
    // The final Newton step after NLDD domain solves re-synchronizes everything
    // via a global BlackoilWellModel::assemble() which performs the full control
    // update pipeline: BlackoilWellModelGeneric::updateAndCommunicateGroupData,
    // BlackoilWellModelGeneric::updateGroupControls,
    // group+individual well constraint checks, network balancing, and gas lift.
    // Importantly, the global step also checks individual constraints (after
    // group constraints), so if a domain solve switched a well to e.g. BHP
    // because of a genuine BHP limit violation, the global step will reach
    // the same conclusion -- the two levels are consistent for individual
    // constraints.

    // Use do_mpi_gather=false to avoid MPI collective operations in domain solves.
    auto loggerGuard = wellModel_.groupStateHelper().pushLogger(/*do_mpi_gather=*/false);

    this->updateWellControls(domain);
    this->assembleWellEq(dt, domain);

    // Update cellRates_ with current contributions from wells in this domain for reservoir linearization
    wellModel_.updateCellRatesForDomain(domain.index, this->well_domain());
}

template<typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
assembleWellEq(const double dt,
               const Domain& domain)
{
    OPM_TIMEBLOCK(assembleWellEq);
    for (const auto& well : wellModel_.localNonshutWells()) {
        if (this->well_domain().at(well->name()) == domain.index) {
            well->assembleWellEq(wellModel_.simulator(),
                                 dt,
                                 wellModel_.groupStateHelper(),
                                 wellModel_.wellState());
        }
    }
}

template<typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
addWellPressureEquations(PressureMatrix& /*jacobian*/,
                         const BVector& /*weights*/,
                         const bool /*use_well_weights*/,
                         const int /*domainIndex*/) const
{
    throw std::logic_error("CPRW is not yet implemented for NLDD subdomains");
    // To fix this function, rdofs should be the size of the domain, and the nw should be the number of wells in the domain
    // int nw = this->numLocalWellsEnd(); // should number of wells in the domain
    // int rdofs = local_num_cells_;  // should be the size of the domain
    // for ( int i = 0; i < nw; i++ ) {
    //     int wdof = rdofs + i;
    //     jacobian[wdof][wdof] = 1.0;// better scaling ?
    // }

    // for ( const auto& well : well_container_ ) {
    //     if (well_domain_.at(well->name()) == domainIndex) {
               // weights should be the size of the domain
    //         well->addWellPressureEquations(jacobian, weights, pressureVarIndex, use_well_weights, this->wellState());
    //     }
    // }
}

template<typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
recoverWellSolutionAndUpdateWellState(const BVector& x,
                                      const int domainIdx)
{
    // Note: no point in trying to do a parallel gathering
    // try/catch here, as this function is not called in
    // parallel but for each individual domain of each rank.
    // Use do_mpi_gather=false to avoid MPI collective operations.
    auto loggerGuard = wellModel_.groupStateHelper().pushLogger(/*do_mpi_gather=*/false);
    for (const auto& well : wellModel_.localNonshutWells()) {
        if (this->well_domain().at(well->name()) == domainIdx) {
            const auto& cells = well->cells();
            x_local_.resize(cells.size());

            for (size_t i = 0; i < cells.size(); ++i) {
                x_local_[i] = x[cells[i]];
            }
            well->recoverWellSolutionAndUpdateWellState(wellModel_.simulator(),
                                                        x_local_,
                                                        wellModel_.groupStateHelper(),
                                                        wellModel_.wellState());
        }
    }
}

template<typename TypeTag>
ConvergenceReport
BlackoilWellModelNldd<TypeTag>::
getWellConvergence(const Domain& domain,
                   const std::vector<Scalar>& B_avg,
                   DeferredLogger& local_deferredLogger) const
{
    const auto& iterCtx = wellModel_.simulator().problem().iterationContext();
    const bool relax_tolerance = iterCtx.shouldRelax(wellModel_.numStrictIterations());

    ConvergenceReport report;
    {
        // Use do_mpi_gather=false to avoid MPI collective operations in domain solves.
        auto loggerGuard = wellModel_.groupStateHelper().pushLogger(/*do_mpi_gather=*/false);

        for (const auto& well : wellModel_.localNonshutWells()) {
            if ((this->well_domain().at(well->name()) == domain.index)) {
                if (well->isOperableAndSolvable() || well->wellIsStopped()) {
                    report += well->getWellConvergence(wellModel_.groupStateHelper(),
                                                       B_avg,
                                                       relax_tolerance);
                } else {
                    ConvergenceReport xreport;
                    using CR = ConvergenceReport;
                    xreport.setWellFailed({CR::WellFailure::Type::Unsolvable,
                                           CR::Severity::Normal, -1, well->name()});
                    report += xreport;
                }
            }
        }
    } // loggerGuard goes out of scope here, before the OpmLog::debug() calls below

    // Log debug messages for NaN or too large residuals.
    if (wellModel_.terminalOutput()) {
        for (const auto& f : report.wellFailures()) {
            if (f.severity() == ConvergenceReport::Severity::NotANumber) {
                local_deferredLogger.debug("NaN residual found with phase " +
                                           std::to_string(f.phase()) +
                                           " for well " + f.wellName());
            } else if (f.severity() == ConvergenceReport::Severity::TooLarge) {
                local_deferredLogger.debug("Too large residual found with phase " +
                                           std::to_string(f.phase()) +
                                           " for well " + f.wellName());
            }
        }
    }
    return report;
}

template<typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
updateWellControls(const Domain& domain)
{
    OPM_TIMEBLOCK(updateWellControls);
    if (!wellModel_.wellsActive()) {
        return;
    }

    // Only individual well constraints (BHP, THP, rate limits) are checked.
    // Group constraints, network balancing, and NUPCOL/VREP updates require
    // MPI collectives and consistent global group state, which are not
    // available during domain-local solves. The final global Newton step
    // handles all group-level control decisions.
    //
    // Individual constraints are physical safety limits that are essential
    // for well convergence (e.g., BHP minimum prevents infeasible operating
    // points). Skipping them would cause domain convergence failures.
    for (const auto& well : wellModel_.localNonshutWells()) {
        if (this->well_domain().at(well->name()) == domain.index) {
            constexpr auto mode = WellInterface<TypeTag>::IndividualOrGroup::Individual;
            well->updateWellControl(wellModel_.simulator(),
                                    mode,
                                    wellModel_.groupStateHelper(),
                                    wellModel_.wellState());
        }
    }
}

template <typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
setupDomains(const std::vector<Domain>& domains)
{
    std::vector<const SubDomainIndices*> genDomains;
    std::ranges::transform(domains, std::back_inserter(genDomains),
                           [](const auto& domain)
                           { return static_cast<const SubDomainIndices*>(&domain); });
    this->calcDomains(genDomains);
}

} // namespace Opm

#endif
