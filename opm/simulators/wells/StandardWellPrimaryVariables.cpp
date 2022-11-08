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
#include <opm/simulators/wells/StandardWellPrimaryVariables.hpp>

#include <dune/common/dynvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

namespace Opm {

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
init()
{
    for (int eqIdx = 0; eqIdx < numWellEq_; ++eqIdx) {
        evaluation_[eqIdx] =
            EvalWell::createVariable(numWellEq_ + Indices::numEq,
                                     value_[eqIdx],
                                     Indices::numEq + eqIdx);

    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
resize(const int numWellEq)
{
    value_.resize(numWellEq, 0.0);
    evaluation_.resize(numWellEq, EvalWell{numWellEq + Indices::numEq, 0.0});
    numWellEq_ = numWellEq;
}


template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
updatePolyMW(const BVectorWell& dwells)
{
    if (well_.isInjector()) {
        for (int perf = 0; perf < well_.numPerfs(); ++perf) {
            const int wat_vel_index = Bhp + 1 + perf;
            const int pskin_index = Bhp + 1 + well_.numPerfs() + perf;

            const double relaxation_factor = 0.9;
            const double dx_wat_vel = dwells[0][wat_vel_index];
            value_[wat_vel_index] -= relaxation_factor * dx_wat_vel;

            const double dx_pskin = dwells[0][pskin_index];
            value_[pskin_index] -= relaxation_factor * dx_pskin;
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
copyToWellStatePolyMW(WellState& well_state) const
{
    if (well_.isInjector()) {
        auto& ws = well_state.well(well_.indexOfWell());
        auto& perf_data = ws.perf_data;
        auto& perf_water_velocity = perf_data.water_velocity;
        auto& perf_skin_pressure = perf_data.skin_pressure;
        for (int perf = 0; perf < well_.numPerfs(); ++perf) {
            perf_water_velocity[perf] = value_[Bhp + 1 + perf];
            perf_skin_pressure[perf] = value_[Bhp + 1 + well_.numPerfs() + perf];
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
copyToWellState(WellState& well_state,
                DeferredLogger& deferred_logger) const
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const PhaseUsage& pu = well_.phaseUsage();
    std::vector<double> F(well_.numPhases(), 0.0);
    [[maybe_unused]] double F_solvent = 0.0;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const int oil_pos = pu.phase_pos[Oil];
        F[oil_pos] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const int water_pos = pu.phase_pos[Water];
            F[water_pos] = value_[WFrac];
            F[oil_pos] -= F[water_pos];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int gas_pos = pu.phase_pos[Gas];
            F[gas_pos] = value_[GFrac];
            F[oil_pos] -= F[gas_pos];
        }

        if constexpr (Indices::enableSolvent) {
            F_solvent = value_[SFrac];
            F[oil_pos] -= F_solvent;
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const int water_pos = pu.phase_pos[Water];
        F[water_pos] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int gas_pos = pu.phase_pos[Gas];
            F[gas_pos] = value_[GFrac];
            F[water_pos] -= F[gas_pos];
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const int gas_pos = pu.phase_pos[Gas];
        F[gas_pos] = 1.0;
    }

    // convert the fractions to be Q_p / G_total to calculate the phase rates
    for (int p = 0; p < well_.numPhases(); ++p) {
        const double scal = well_.scalingFactor(p);
        // for injection wells, there should only one non-zero scaling factor
        if (scal > 0) {
            F[p] /= scal ;
        } else {
            // this should only happens to injection wells
            F[p] = 0.;
        }
    }

    // F_solvent is added to F_gas. This means that well_rate[Gas] also contains solvent.
    // More testing is needed to make sure this is correct for well groups and THP.
    if constexpr (Indices::enableSolvent) {
        F_solvent /= well_.scalingFactor(Indices::contiSolventEqIdx);
        F[pu.phase_pos[Gas]] += F_solvent;
    }

    auto& ws = well_state.well(well_.indexOfWell());
    ws.bhp = value_[Bhp];

    // calculate the phase rates based on the primary variables
    // for producers, this is not a problem, while not sure for injectors here
    if (well_.isProducer()) {
        const double g_total = value_[WQTotal];
        for (int p = 0; p < well_.numPhases(); ++p) {
            ws.surface_rates[p] = g_total * F[p];
        }
    } else { // injectors
        for (int p = 0; p < well_.numPhases(); ++p) {
            ws.surface_rates[p] = 0.0;
        }
        switch (well_.wellEcl().injectorType()) {
        case InjectorType::WATER:
            ws.surface_rates[pu.phase_pos[Water]] = value_[WQTotal];
            break;
        case InjectorType::GAS:
            ws.surface_rates[pu.phase_pos[Gas]] = value_[WQTotal];
            break;
        case InjectorType::OIL:
            ws.surface_rates[pu.phase_pos[Oil]] = value_[WQTotal];
            break;
        case InjectorType::MULTI:
            // Not supported.
            deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                                    "Multi phase injectors are not supported, requested for well " + well_.name());
            break;
        }
    }
}

#define INSTANCE(...) \
template class StandardWellPrimaryVariables<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)

// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)

}
