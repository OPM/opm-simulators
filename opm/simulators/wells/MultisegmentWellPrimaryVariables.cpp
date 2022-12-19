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
#include <opm/simulators/wells/MultisegmentWellPrimaryVariables.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <algorithm>

namespace Opm {

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
resize(const int numSegments)
{
    value_.resize(numSegments);
    evaluation_.resize(numSegments);
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
init()
{
    for (size_t seg = 0; seg < value_.size(); ++seg) {
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            evaluation_[seg][eq_idx] = 0.0;
            evaluation_[seg][eq_idx].setValue(value_[seg][eq_idx]);
            evaluation_[seg][eq_idx].setDerivative(eq_idx + Indices::numEq, 1.0);
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
update(const WellState& well_state)
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Gas = BlackoilPhases::Vapour;

    // TODO: to test using rate conversion coefficients to see if it will be better than
    // this default one
    if (!well_.isOperableAndSolvable() && !well_.wellIsStopped())
        return;

    const Well& well = well_.wellEcl();

    // the index of the top segment in the WellState
    const auto& ws = well_state.well(well_.indexOfWell());
    const auto& segments = ws.segments;
    // maybe a throw for parallel running?
    assert(segments.size() == value_.size());
    const auto& segment_rates = segments.rates;
    const auto& segment_pressure = segments.pressure;
    const PhaseUsage& pu = well_.phaseUsage();

    for (size_t seg = 0; seg < value_.size(); ++seg) {
        // calculate the total rate for each segment
        double total_seg_rate = 0.0;
        // the segment pressure
        value_[seg][SPres] = segment_pressure[seg];
        // TODO: under what kind of circustances, the following will be wrong?
        // the definition of g makes the gas phase is always the last phase
        for (int p = 0; p < well_.numPhases(); p++) {
            total_seg_rate += well_.scalingFactor(p) * segment_rates[well_.numPhases() * seg + p];
        }

        if (seg == 0) {
            if (well_.isInjector()) {
                total_seg_rate = std::max(total_seg_rate, 0.);
            } else {
                total_seg_rate = std::min(total_seg_rate, 0.);
            }
        }
        value_[seg][WQTotal] = total_seg_rate;
        if (std::abs(total_seg_rate) > 0.) {
            if (has_wfrac_variable) {
                const int water_pos = pu.phase_pos[Water];
                value_[seg][WFrac] = well_.scalingFactor(water_pos) * segment_rates[well_.numPhases() * seg + water_pos] / total_seg_rate;
            }
            if (has_gfrac_variable) {
                const int gas_pos = pu.phase_pos[Gas];
                value_[seg][GFrac] = well_.scalingFactor(gas_pos) * segment_rates[well_.numPhases() * seg + gas_pos] / total_seg_rate;
            }
        } else { // total_seg_rate == 0
            if (well_.isInjector()) {
                // only single phase injection handled
                auto phase = well.getInjectionProperties().injectorType;

                if (has_wfrac_variable) {
                    if (phase == InjectorType::WATER) {
                        value_[seg][WFrac] = 1.0;
                    } else {
                        value_[seg][WFrac] = 0.0;
                    }
                }

                if (has_gfrac_variable) {
                    if (phase == InjectorType::GAS) {
                        value_[seg][GFrac] = 1.0;
                    } else {
                        value_[seg][GFrac] = 0.0;
                    }
                }

            } else if (well_.isProducer()) { // producers
                if (has_wfrac_variable) {
                    value_[seg][WFrac] = 1.0 / well_.numPhases();
                }

                if (has_gfrac_variable) {
                    value_[seg][GFrac] = 1.0 / well_.numPhases();
                }
            }
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
updateNewton(const BVectorWell& dwells,
             const double relaxation_factor,
             const double dFLimit,
             const double max_pressure_change)
{
    const std::vector<std::array<double, numWellEq>> old_primary_variables = value_;

    for (size_t seg = 0; seg < value_.size(); ++seg) {
        if (has_wfrac_variable) {
            const int sign = dwells[seg][WFrac] > 0. ? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][WFrac]) * relaxation_factor, dFLimit);
            value_[seg][WFrac] = old_primary_variables[seg][WFrac] - dx_limited;
        }

        if (has_gfrac_variable) {
            const int sign = dwells[seg][GFrac] > 0. ? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][GFrac]) * relaxation_factor, dFLimit);
            value_[seg][GFrac] = old_primary_variables[seg][GFrac] - dx_limited;
        }

        // handling the overshooting or undershooting of the fractions
        this->processFractions(seg);

        // update the segment pressure
        {
            const int sign = dwells[seg][SPres] > 0.? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][SPres]) * relaxation_factor, max_pressure_change);
            value_[seg][SPres] = std::max(old_primary_variables[seg][SPres] - dx_limited, 1e5);
        }

        // update the total rate // TODO: should we have a limitation of the total rate change?
        {
            value_[seg][WQTotal] = old_primary_variables[seg][WQTotal] - relaxation_factor * dwells[seg][WQTotal];

            // make sure that no injector produce and no producer inject
            if (seg == 0) {
                if (well_.isInjector()) {
                    value_[seg][WQTotal] = std::max(value_[seg][WQTotal], 0.0);
                } else {
                    value_[seg][WQTotal] = std::min(value_[seg][WQTotal], 0.0);
                }
            }
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
processFractions(const int seg)
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const PhaseUsage& pu = well_.phaseUsage();

    std::vector<double> fractions(well_.numPhases(), 0.0);

    assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
    const int oil_pos = pu.phase_pos[Oil];
    fractions[oil_pos] = 1.0;

    if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
        const int water_pos = pu.phase_pos[Water];
        fractions[water_pos] = value_[seg][WFrac];
        fractions[oil_pos] -= fractions[water_pos];
    }

    if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
        const int gas_pos = pu.phase_pos[Gas];
        fractions[gas_pos] = value_[seg][GFrac];
        fractions[oil_pos] -= fractions[gas_pos];
    }

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const int water_pos = pu.phase_pos[Water];
        if (fractions[water_pos] < 0.0) {
            if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
                fractions[pu.phase_pos[Gas]] /= (1.0 - fractions[water_pos]);
            }
            fractions[oil_pos] /= (1.0 - fractions[water_pos]);
            fractions[water_pos] = 0.0;
        }
    }

    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const int gas_pos = pu.phase_pos[Gas];
        if (fractions[gas_pos] < 0.0) {
            if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
                fractions[pu.phase_pos[Water]] /= (1.0 - fractions[gas_pos]);
            }
            fractions[oil_pos] /= (1.0 - fractions[gas_pos]);
            fractions[gas_pos] = 0.0;
        }
    }

    if (fractions[oil_pos] < 0.0) {
        if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
            fractions[pu.phase_pos[Water]] /= (1.0 - fractions[oil_pos]);
        }
        if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            fractions[pu.phase_pos[Gas]] /= (1.0 - fractions[oil_pos]);
        }
        fractions[oil_pos] = 0.0;
    }

    if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
        value_[seg][WFrac] = fractions[pu.phase_pos[Water]];
    }

    if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
        value_[seg][GFrac] = fractions[pu.phase_pos[Gas]];
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
volumeFraction(const int seg,
               const unsigned compIdx) const
{
    if (has_wfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)) {
        return evaluation_[seg][WFrac];
    }

    if (has_gfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
        return evaluation_[seg][GFrac];
    }

    // Oil fraction
    EvalWell oil_fraction = 1.0;
    if (has_wfrac_variable) {
        oil_fraction -= evaluation_[seg][WFrac];
    }

    if (has_gfrac_variable) {
        oil_fraction -= evaluation_[seg][GFrac];
    }
    /* if (has_solvent) {
        oil_fraction -= evaluation_[seg][SFrac];
    } */
    return oil_fraction;
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
volumeFractionScaled(const int seg,
                     const int comp_idx) const
{
    // For reservoir rate control, the distr in well control is used for the
    // rate conversion coefficients. For the injection well, only the distr of the injection
    // phase is not zero.
    const double scale = well_.scalingFactor(well_.ebosCompIdxToFlowCompIdx(comp_idx));
    if (scale > 0.) {
        return this->volumeFraction(seg, comp_idx) / scale;
    }

    return this->volumeFraction(seg, comp_idx);
}

#define INSTANCE(...) \
template class MultisegmentWellPrimaryVariables<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
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
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,1u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,2u,0u>)

}
