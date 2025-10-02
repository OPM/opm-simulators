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
template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
template<typename Scalar, typename IndexTraits> class WellState;

template<typename Scalar, typename IndexTraits>
class WellConvergence
{
public:
    explicit WellConvergence(const WellInterfaceGeneric<Scalar, IndexTraits>& well)
        : well_(well)
    {}

    struct Tolerances {
        Scalar bhp; //!< Tolerance for bhp controlled well
        Scalar thp; //!< Tolerance for thp controlled well
        Scalar rates; //!< Tolerance for a rate controlled well
        Scalar grup; //!< Tolerance for a grup controlled well
        Scalar max_residual_allowed; //!< Max residual allowd
    };

    // checking the convergence of the well control equations
    void checkConvergenceControlEq(const WellState<Scalar, IndexTraits>& well_state,
                                   const Tolerances& tolerances,
                                   const Scalar well_control_residual,
                                   const bool well_is_stopped,
                                   ConvergenceReport& report,
                                   DeferredLogger& deferred_logger) const;

    void checkConvergencePolyMW(const std::vector<Scalar>& res,
                                const int Bhp,
                                const Scalar maxResidualAllowed,
                                ConvergenceReport& report) const;

private:
    const WellInterfaceGeneric<Scalar, IndexTraits>& well_;
};

}

#endif // OPM_WELL_CONVERGENCE_HEADER_INCLUDED
