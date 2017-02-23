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

#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <vector>

namespace Opm {

    // Forward declarations
    class SimulationDataContainer;
    class WellStateFullyImplicitBlackoil;
    class WellStateFullyImplicitBlackoilDense;

    std::vector< double > destripe( const std::vector< double >& v,
                                    size_t stride,
                                    size_t offset );

    std::vector< double >& stripe( const std::vector< double >& v,
                                   size_t stride,
                                   size_t offset,
                                   std::vector< double >& dst );

    data::Solution simToSolution( const SimulationDataContainer& reservoir,
                                  PhaseUsage phases );

    void solutionToSim( const data::Solution& sol,
                        PhaseUsage phases,
                        SimulationDataContainer& state );

    void wellsToState( const data::Wells& wells,
                       PhaseUsage phases,
                       WellStateFullyImplicitBlackoil& state );

    void wellsToState( const data::Wells& wells,
                       PhaseUsage phases,
                       WellStateFullyImplicitBlackoilDense& state );

}

#endif //OPM_SIMULATORS_COMPAT_HPP
