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


#ifndef OPM_WELL_CONSTRAINTS_HEADER_INCLUDED
#define OPM_WELL_CONSTRAINTS_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <functional>
#include <vector>
#include <optional>

namespace Opm
{

class DeferredLogger;
using RegionId = int;
class Rates;
template<typename Scalar, typename IndexTraits> class SingleWellState;
class SummaryState;
template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
enum class WellInjectorCMode;
enum class WellProducerCMode;

//! \brief Class for computing well group constraints.
template<typename Scalar, typename IndexTraits>
class WellConstraints {
public:
    //! \brief Constructor sets reference to well.
    explicit WellConstraints(const WellInterfaceGeneric<Scalar, IndexTraits>& well) : well_(well) {}

    using RateConvFunc = std::function<void(const RegionId, const int,
                                            const std::vector<Scalar>&,
                                            std::vector<Scalar>&)>;

    bool
    checkIndividualConstraints(SingleWellState<Scalar, IndexTraits>& ws,
                               const SummaryState& summaryState,
                               const RateConvFunc& calcReservoirVoidageRates,
                               bool& thp_limit_violated_but_not_switched,
                               DeferredLogger& deferred_logger,
                               const std::optional<Well::InjectionControls>& inj_controls = std::nullopt,
                               const std::optional<Well::ProductionControls>& prod_controls = std::nullopt) const;

private:
    WellInjectorCMode
    activeInjectionConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                              const SummaryState& summaryState,
                              bool& thp_limit_violated_but_not_switched,
                              DeferredLogger& deferred_logger,
                              const std::optional<Well::InjectionControls>& inj_controls = std::nullopt) const;

    WellProducerCMode
    activeProductionConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                               const SummaryState& summaryState,
                               const RateConvFunc& calcReservoirVoidageRates,
                               bool& thp_limit_violated_but_not_switched,
                               DeferredLogger& deferred_logger,
                               const std::optional<Well::ProductionControls>& prod_controls = std::nullopt) const;

    const WellInterfaceGeneric<Scalar, IndexTraits>& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_CONSTRAINTS_HEADER_INCLUDED
