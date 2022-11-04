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


#ifndef OPM_STANDARDWELL_GENERIC_HEADER_INCLUDED
#define OPM_STANDARDWELL_GENERIC_HEADER_INCLUDED

#include <opm/simulators/wells/WellHelpers.hpp>

#include <optional>
#include <vector>
#include <functional>

namespace Opm
{

class ConvergenceReport;
class DeferredLogger;
class ParallelWellInfo;
class Schedule;
class SummaryState;
class WellInterfaceGeneric;
class WellState;

template<class Scalar>
class StandardWellGeneric
{
protected:
    StandardWellGeneric(const WellInterfaceGeneric& baseif);

    // calculate a relaxation factor to avoid overshoot of total rates
    static double relaxationFactorRate(const std::vector<double>& primary_variables,
                                       const double newton_update);

    // relaxation factor considering only one fraction value
    static double relaxationFactorFraction(const double old_value,
                                           const double dx);

    void computeConnectionPressureDelta();

    // Base interface reference
    const WellInterfaceGeneric& baseif_;

    // densities of the fluid in each perforation
    std::vector<double> perf_densities_;
    // pressure drop between different perforations
    std::vector<double> perf_pressure_diffs_;

    double getRho() const
    {
        return this->perf_densities_.empty() ? 0.0 : perf_densities_[0];
    }
};

}

#endif // OPM_STANDARDWELL_GENERIC_HEADER_INCLUDED
