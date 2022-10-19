/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELL_BPH_THP_CALCULATOR_HEADER_INCLUDED
#define OPM_WELL_BPH_THP_CALCULATOR_HEADER_INCLUDED

#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace Opm
{

class DeferredLogger;
class SummaryState;
class WellInterfaceGeneric;

//! \brief Class for computing BHP limits.
class WellBhpThpCalculator {
public:
    //! \brief Constructor sets reference to well.
    WellBhpThpCalculator(const WellInterfaceGeneric& well) : well_(well) {}

    //! \brief Checks if well has THP constraints.
    bool wellHasTHPConstraints(const SummaryState& summaryState) const;

    //! \brief Get THP constraint for well.
    double getTHPConstraint(const SummaryState& summaryState) const;

    //! \brief Obtain the most strict BHP from BHP limits.
    double mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const;

    //! \brief Calculates THP from BHP.
    double calculateThpFromBhp(const std::vector<double>& rates,
                               const double bhp,
                               const double rho,
                               const double alq,
                               DeferredLogger& deferred_logger) const;

    //! \brief Compute BHP from THP limit for a producer.
    std::optional<double>
    computeBhpAtThpLimitProd(const std::function<std::vector<double>(const double)>& frates,
                             const SummaryState& summary_state,
                             const double maxPerfPress,
                             const double rho,
                             const double alq_value,
                             DeferredLogger& deferred_logger) const;

private:
    //! \brief Calculate max BHP.
    std::optional<double>
    bhpMax(const std::function<double(const double)>& fflo,
           const double bhp_limit,
           const double maxPerfPress,
           const double vfp_flo_front,
           DeferredLogger& deferred_logger) const;

    //! \brief Common code for finding BHP from THP limit for producers/injectors.
    std::optional<double>
    computeBhpAtThpLimit(const std::function<std::vector<double>(const double)>& frates,
                         const std::function<double(const std::vector<double>)>& fbhp,
                         const std::array<double, 2>& range,
                         DeferredLogger& deferred_logger) const;

    //! \brief Find limits using bisection.
    bool bisectBracket(const std::function<double(const double)>& eq,
                       const std::array<double, 2>& range,
                       double& low, double& high,
                       std::optional<double>& approximate_solution,
                       DeferredLogger& deferred_logger) const;

    //! \brief Find limits using brute-force solver.
    static bool bruteForceBracket(const std::function<double(const double)>& eq,
                                  const std::array<double, 2>& range,
                                  double& low, double& high,
                                  DeferredLogger& deferred_logger);

    const WellInterfaceGeneric& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_BHP_THP_CALCULATOR_HEADER_INCLUDED
