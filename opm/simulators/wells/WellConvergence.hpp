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


#ifndef OPM_WELL_CONVERGENCE_HEADER_INCLUDED
#define OPM_WELL_CONVERGENCE_HEADER_INCLUDED

#include <vector>

namespace Opm
{

class ConvergenceReport;
class DeferredLogger;
class WellInterfaceGeneric;
class WellState;

class WellConvergence
{
public:
    WellConvergence(const WellInterfaceGeneric& well)
        : well_(well)
    {}

    struct Tolerances {
        double bhp; //!< Tolerance for bhp controlled well
        double thp; //!< Tolerance for thp controlled well
        double rates; //!< Tolerance for a rate controlled well
        double grup; //!< Tolerance for a grup controlled well
        double max_residual_allowed; //!< Max residual allowd
    };

    // checking the convergence of the well control equations
    void checkConvergenceControlEq(const WellState& well_state,
                                   const Tolerances& tolerances,
                                   const double well_control_residual,
                                   const bool well_is_stopped, 
                                   ConvergenceReport& report,
                                   DeferredLogger& deferred_logger) const;

    void checkConvergencePolyMW(const std::vector<double>& res,
                                const int Bhp,
                                const double maxResidualAllowed,
                                ConvergenceReport& report) const;

private:
    const WellInterfaceGeneric& well_;
};

}

#endif // OPM_WELL_CONVERGENCE_HEADER_INCLUDED
