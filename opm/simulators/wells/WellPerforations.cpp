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
#include <opm/simulators/wells/WellPerforations.hpp>

#include <opm/material/densead/Evaluation.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/wells/PerforationData.hpp>

namespace Opm {

template<class FluidSystem, class Indices, class Scalar, class Value>
void WellPerforations<FluidSystem,Indices,Scalar,Value>::
gasOilRateProd(std::vector<Value>& cq_s,
               PerforationRates& perf_rates,
               const Value& rv,
               const Value& rs,
               const Value& rvw) const
{
    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    const Value cq_sOil = cq_s[oilCompIdx];
    const Value cq_sGas = cq_s[gasCompIdx];
    const Value dis_gas = rs * cq_sOil;
    const Value vap_oil = rv * cq_sGas;

    cq_s[gasCompIdx] += dis_gas;
    cq_s[oilCompIdx] += vap_oil;

    // recording the perforation solution gas rate and solution oil rates
    if (well_.isProducer()) {
        perf_rates.dis_gas = getValue(dis_gas);
        perf_rates.vap_oil = getValue(vap_oil);
    }

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        const Value vap_wat = rvw * cq_sGas;
        cq_s[waterCompIdx] += vap_wat;
        if (well_.isProducer()) {
            perf_rates.vap_wat = getValue(vap_wat);
        }
    }
}

#define INSTANCE(Dim, ...) \
template class WellPerforations<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>, \
                                __VA_ARGS__, double, double>; \
template class WellPerforations<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>, \
                                __VA_ARGS__, double, \
                                DenseAd::Evaluation<double,-1,Dim>>;

INSTANCE(4u, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(5u, BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(9u, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(8u, BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)

INSTANCE(8u, BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(10u, BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(10u, BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)

} // namespace Opm
