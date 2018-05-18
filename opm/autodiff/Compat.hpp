/*
   Copyright 2016 Statoil ASA
   2016 IRIS

   This file is part of the Open Porous Media project (OPM).

   OPM is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OPM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with OPM. If not, see <http://www.gnu.org/licenses/>.
   */

#ifndef OPM_SIMULATORS_COMPAT_HPP
#define OPM_SIMULATORS_COMPAT_HPP

#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <vector>

namespace Opm {

    // Forward declarations
    class SimulationDataContainer;
    class WellStateFullyImplicitBlackoil;

    /// Extract single data vector from striped data.
    /// \return   u such that u[i] = v[offset + i*stride].
    std::vector< double > destripe( const std::vector< double >& v,
                                    size_t stride,
                                    size_t offset );

    /// Inject single data vector into striped data.
    /// \return   reference to dst input, that is changed so that
    ///           dst[offset + i*stride] = v[i]. This is done for
    ///           i = 0..(dst.size()/stride).
    std::vector< double >& stripe( const std::vector< double >& v,
                                   size_t stride,
                                   size_t offset,
                                   std::vector< double >& dst );

    /// Returns Solution with the following fields:
    ///   PRESSURE, TEMP (unconditionally)
    ///   SWAT, SGAS, RS, RV, SSOL (if appropriate fields present in input)
    /// If use_si_units is true, the fields will have the measure UnitSystem::measure::identity,
    /// and therefore *not* be converted to customary units (depending on family) upon output.
    data::Solution simToSolution( const SimulationDataContainer& reservoir,
                                  const bool use_si_units,
                                  PhaseUsage phases );

    /// Copies the following fields from sol into state (all conditionally):
    ///   PRESSURE, TEMP, SWAT, SGAS, RS, RV, SSOL
    /// Also handles extra data such as hysteresis parameters, SOMAX, etc.
    void solutionToSim( const RestartValue& restart_value,
                        PhaseUsage phases,
                        SimulationDataContainer& state );

    /// Copies the following fields from wells into state.
    ///   bhp, temperature, currentControls, wellRates, perfPress, perfRates, perfPhaseRates
    void wellsToState( const data::Wells& wells,
                       PhaseUsage phases,
                       WellStateFullyImplicitBlackoil& state );

}

#endif //OPM_SIMULATORS_COMPAT_HPP
