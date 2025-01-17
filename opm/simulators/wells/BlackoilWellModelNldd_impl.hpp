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
assemble(const int /*iterationIdx*/,
         const double dt,
         const Domain& domain)
{
    // We assume that calculateExplicitQuantities() and
    // prepareTimeStep() have been called already for the entire
    // well model, so we do not need to do it here (when
    // iterationIdx is 0).

    DeferredLogger local_deferredLogger;
    this->updateWellControls(local_deferredLogger, domain);
    this->initPrimaryVariablesEvaluation(domain);
    this->assembleWellEq(dt, domain, local_deferredLogger);
}

template<typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
assembleWellEq(const double dt,
               const Domain& domain,
               DeferredLogger& deferred_logger)
{
    for (const auto& well : wellModel_.localNonshutWells()) {
        if (this->well_domain().at(well->name()) == domain.index) {
            well->assembleWellEq(wellModel_.simulator(),
                                 dt,
                                 wellModel_.wellState(),
                                 wellModel_.groupState(),
                                 deferred_logger);
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
    DeferredLogger local_deferredLogger;
    for (const auto& well : wellModel_.localNonshutWells()) {
        if (this->well_domain().at(well->name()) == domainIdx) {
            const auto& cells = well->cells();
            x_local_.resize(cells.size());

            for (size_t i = 0; i < cells.size(); ++i) {
                x_local_[i] = x[cells[i]];
            }
            well->recoverWellSolutionAndUpdateWellState(wellModel_.simulator(),
                                                        x_local_,
                                                        wellModel_.wellState(),
                                                        local_deferredLogger);
        }
    }
    // TODO: avoid losing the logging information that could
    // be stored in the local_deferredlogger in a parallel case.
    if (wellModel_.terminalOutput()) {
        local_deferredLogger.logMessages();
    }
}

template<typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
initPrimaryVariablesEvaluation(const Domain& domain) const
{
    for (auto& well : wellModel_.localNonshutWells()) {
        if (this->well_domain().at(well->name()) == domain.index) {
            well->initPrimaryVariablesEvaluation();
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
    const int iterationIdx = wellModel_.simulator().model().newtonMethod().numIterations();
    const bool relax_tolerance = iterationIdx > wellModel_.numStrictIterations();

    ConvergenceReport report;
    for (const auto& well : wellModel_.localNonshutWells()) {
        if ((this->well_domain().at(well->name()) == domain.index)) {
            if (well->isOperableAndSolvable() || well->wellIsStopped()) {
                report += well->getWellConvergence(wellModel_.simulator(),
                                                   wellModel_.wellState(),
                                                   B_avg,
                                                   local_deferredLogger,
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
updateWellControls(DeferredLogger& deferred_logger,
                   const Domain& domain)
{
    if (!wellModel_.wellsActive()) {
        return;
    }

    // TODO: decide on and implement an approach to handling of
    // group controls, network and similar for domain solves.

    // Check only individual well constraints and communicate.
    for (const auto& well : wellModel_.localNonshutWells()) {
        if (this->well_domain().at(well->name()) == domain.index) {
            constexpr auto mode = WellInterface<TypeTag>::IndividualOrGroup::Individual;
            well->updateWellControl(wellModel_.simulator(),
                                    mode,
                                    wellModel_.wellState(),
                                    wellModel_.groupState(),
                                    deferred_logger);
        }
    }
}

template <typename TypeTag>
void
BlackoilWellModelNldd<TypeTag>::
setupDomains(const std::vector<Domain>& domains)
{
    std::vector<const SubDomainIndices*> genDomains;
    std::transform(domains.begin(), domains.end(),
                   std::back_inserter(genDomains),
                   [](const auto& domain)
                   { return static_cast<const SubDomainIndices*>(&domain); });
    this->calcDomains(genDomains);
}

} // namespace Opm

#endif
