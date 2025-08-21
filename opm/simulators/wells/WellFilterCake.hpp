/*
  Copyright 2023 Equinor

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

#ifndef OPM_WELL_FILTER_CAKE_HEADER_INCLUDED
#define OPM_WELL_FILTER_CAKE_HEADER_INCLUDED

#include <cstddef>
#include <vector>

namespace Opm {

class DeferredLogger;
template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
template<typename Scalar, typename IndexTraits> class WellState;

//! \brief Class for well calculations related to filter cakes.
template<typename Scalar, typename IndexTraits>
class WellFilterCake {
public:
    //! \brief Post-step filtration model updates
    //! \details Calculates the filtrate deposition volumes and associated skin factors / injectivity multipliers
    void updatePostStep(const WellInterfaceGeneric<Scalar, IndexTraits>& well,
                WellState<Scalar, IndexTraits>& well_state,
                const double dt,
                const Scalar conc,
                const std::size_t water_index,
                DeferredLogger& deferred_logger);

    //! \brief Pre-step filtration model updates
    //! \details Applies filter cake cleaning
    void updatePreStep(const WellInterfaceGeneric<Scalar, IndexTraits>& well,
                       DeferredLogger& deferred_logger);

    //! \brief Returns a const-ref to multipliers.
    const std::vector<Scalar>& multipliers() const { return inj_fc_multiplier_; }

private:
    //! \brief Update the multiplier for well transmissbility due to cake filtration.
    void updateSkinFactorsAndMultipliers(const WellInterfaceGeneric<Scalar, IndexTraits>& well,
                                         WellState<Scalar, IndexTraits>& well_state,
                                         const double dt,
                                         const std::size_t water_index,
                                         DeferredLogger& deferred_logger);
    template<class Conn>
    void updateMultiplier(const Conn& conn, const int perf);

    //! \brief Apply cleaning multipliers to skin factors and reduce cake thickness accordingly
    //! \details The cake thickness is re-computed to give the new (reduced) skin factor with current cake properties
    void applyCleaning(const WellInterfaceGeneric<Scalar, IndexTraits>& well,
                       DeferredLogger& deferred_logger);


    std::vector<Scalar> inj_fc_multiplier_; //!< Multiplier due to injection filtration cake
    std::vector<Scalar> skin_factor_;
    std::vector<Scalar> thickness_;
};

}

#endif // OPM_WELL_FILTER_CAKE_HEADER_INCLUDED
