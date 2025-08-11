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

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/wells/WellConvergence.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cmath>
#include <stdexcept>

namespace Opm {

template<typename FluidSystem, typename Indices>
void WellConvergence<FluidSystem, Indices>::
checkConvergenceControlEq(const WellState<FluidSystem, Indices>& well_state,
                          const Tolerances& tolerances,
                          const Scalar well_control_residual,
                          const bool well_is_stopped, 
                          ConvergenceReport& report,
                          DeferredLogger& deferred_logger) const
{
    Scalar control_tolerance = 0.;
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

template<typename FluidSystem, typename Indices>
void WellConvergence<FluidSystem, Indices>::
checkConvergencePolyMW(const std::vector<Scalar>& res,
                       const int Bhp,
                       const Scalar maxResidualAllowed,
                       ConvergenceReport& report) const
{
  if (well_.isInjector()) {
      //  checking the convergence of the perforation rates
      const Scalar wat_vel_tol = 1.e-8;
      const int dummy_component = -1;
      using CR = ConvergenceReport;
      const auto wat_vel_failure_type = CR::WellFailure::Type::MassBalance;
      for (int perf = 0; perf < well_.numLocalPerfs(); ++perf) {
          const Scalar wat_vel_residual = res[Bhp + 1 + perf];
          if (std::isnan(wat_vel_residual)) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::NotANumber, dummy_component, well_.name()});
          } else if (wat_vel_residual > maxResidualAllowed * 10.) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::TooLarge, dummy_component, well_.name()});
          } else if (wat_vel_residual > wat_vel_tol) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::Normal, dummy_component, well_.name()});
          }
      }

      // checking the convergence of the skin pressure
      const Scalar pskin_tol = 1000.; // 1000 pascal
      const auto pskin_failure_type = CR::WellFailure::Type::Pressure;
      for (int perf = 0; perf < well_.numLocalPerfs(); ++perf) {
          const Scalar pskin_residual = res[Bhp + 1 + perf + well_.numLocalPerfs()];
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

    template<class Scalar>
    using FS = BlackOilFluidSystem<Scalar, BlackOilDefaultFluidSystemIndices>;

#define INSTANTIATE(T,...) \
    template class WellConvergence<FS<T>, __VA_ARGS__>;

#define INSTANTIATE_TYPE(T)                                                  \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,false,1u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,true,false,0u,0u>)             \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,true,0u,0u>)             \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,true,2u,0u>)             \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<1u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,1u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,1u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,1u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,1u,false,true,0u,0u>)             \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<1u,0u,0u,0u,true,false,0u,0u>)

    INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
    INSTANTIATE_TYPE(float)
#endif

}
