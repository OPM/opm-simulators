/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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

#include <config.h>
#include <opm/simulators/wells/WellConvergence.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cmath>
#include <stdexcept>

namespace Opm
{

void WellConvergence::
checkConvergenceControlEq(const WellState& well_state,
                          const Tolerances& tolerances,
                          const double well_control_residual,
                          const bool well_is_stopped, 
                          ConvergenceReport& report,
                          DeferredLogger& deferred_logger) const
{
    double control_tolerance = 0.;
    using CR = ConvergenceReport;
    CR::WellFailure::Type ctrltype = CR::WellFailure::Type::Invalid;

    const int well_index = well_.indexOfWell();
    const auto& ws = well_state.well(well_index);
    if (well_is_stopped) {
        ctrltype = CR::WellFailure::Type::ControlRate;
        control_tolerance = tolerances.rates; // use smaller tolerance for zero control?
    }
    else if (well_.isInjector() )
    {
        auto current = ws.injection_cmode;
        switch(current) {
        case Well::InjectorCMode::THP:
            ctrltype = CR::WellFailure::Type::ControlTHP;
            control_tolerance = tolerances.thp;
            break;
        case Well::InjectorCMode::BHP:
            ctrltype = CR::WellFailure::Type::ControlBHP;
            control_tolerance = tolerances.bhp;
            break;
        case Well::InjectorCMode::RATE:
        case Well::InjectorCMode::RESV:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = tolerances.rates;
            break;
        case Well::InjectorCMode::GRUP:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = tolerances.grup;
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error,
                             "Unknown well control control types for well " + well_.name(),
                             deferred_logger);
        }
    }
    else if (well_.isProducer() )
    {
        auto current = ws.production_cmode;
        switch(current) {
        case Well::ProducerCMode::THP:
            ctrltype = CR::WellFailure::Type::ControlTHP;
            control_tolerance = tolerances.thp;
            break;
        case Well::ProducerCMode::BHP:
            ctrltype = CR::WellFailure::Type::ControlBHP;
            control_tolerance = tolerances.bhp;
            break;
        case Well::ProducerCMode::ORAT:
        case Well::ProducerCMode::WRAT:
        case Well::ProducerCMode::GRAT:
        case Well::ProducerCMode::LRAT:
        case Well::ProducerCMode::RESV:
        case Well::ProducerCMode::CRAT:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = tolerances.rates;
            break;
        case Well::ProducerCMode::GRUP:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = tolerances.grup;
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error,
                             "Unknown well control control types for well " + well_.name(),
                             deferred_logger);
        }
    }

    const int dummy_component = -1;
    if (std::isnan(well_control_residual)) {
        report.setWellFailed({ctrltype, CR::Severity::NotANumber, dummy_component, well_.name()});
    } else if (well_control_residual > tolerances.max_residual_allowed * 10.) {
        report.setWellFailed({ctrltype, CR::Severity::TooLarge, dummy_component, well_.name()});
    } else if (well_control_residual > control_tolerance) {
        report.setWellFailed({ctrltype, CR::Severity::Normal, dummy_component, well_.name()});
    }
}

void
WellConvergence::
checkConvergencePolyMW(const std::vector<double>& res,
                       const int Bhp,
                       const double maxResidualAllowed,
                       ConvergenceReport& report) const
{
  if (well_.isInjector()) {
      //  checking the convergence of the perforation rates
      const double wat_vel_tol = 1.e-8;
      const int dummy_component = -1;
      using CR = ConvergenceReport;
      const auto wat_vel_failure_type = CR::WellFailure::Type::MassBalance;
      for (int perf = 0; perf < well_.numPerfs(); ++perf) {
          const double wat_vel_residual = res[Bhp + 1 + perf];
          if (std::isnan(wat_vel_residual)) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::NotANumber, dummy_component, well_.name()});
          } else if (wat_vel_residual > maxResidualAllowed * 10.) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::TooLarge, dummy_component, well_.name()});
          } else if (wat_vel_residual > wat_vel_tol) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::Normal, dummy_component, well_.name()});
          }
      }

      // checking the convergence of the skin pressure
      const double pskin_tol = 1000.; // 1000 pascal
      const auto pskin_failure_type = CR::WellFailure::Type::Pressure;
      for (int perf = 0; perf < well_.numPerfs(); ++perf) {
          const double pskin_residual = res[Bhp + 1 + perf + well_.numPerfs()];
          if (std::isnan(pskin_residual)) {
              report.setWellFailed({pskin_failure_type, CR::Severity::NotANumber, dummy_component, well_.name()});
          } else if (pskin_residual > maxResidualAllowed * 10.) {
              report.setWellFailed({pskin_failure_type, CR::Severity::TooLarge, dummy_component, well_.name()});
          } else if (pskin_residual > pskin_tol) {
              report.setWellFailed({pskin_failure_type, CR::Severity::Normal, dummy_component, well_.name()});
          }
      }
  }
}

}
