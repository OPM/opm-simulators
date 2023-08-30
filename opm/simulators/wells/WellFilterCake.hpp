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
class WellInterfaceGeneric;
class WellState;

//! \brief Class for well calculations related to filter cakes.
class WellFilterCake {
public:
    //! \brief Update the water injection volume.
    //! \details Used for calculation related to cake filtration due to injection activity.
    void updateFiltrationParticleVolume(const WellInterfaceGeneric& well,
                                        const double dt,
                                        const double conc,
                                        const std::size_t water_index,
                                        WellState& well_state);

    //! \brief Update the multiplier for well transmissbility due to cake filtration.
    void updateInjFCMult(const WellInterfaceGeneric& well,
                         WellState& well_state,
                         DeferredLogger& deferred_logger);

    //! \brief Returns a const-ref to multipliers.
    const std::vector<double>& multipliers() const
    {
        return inj_fc_multiplier_;
    }

private:
    std::vector<double> filtration_particle_volume_; //!<// Volume of filtration particles during water injection
    std::vector<double> inj_fc_multiplier_; //!< Multiplier due to injection filtration cake
};

}

#endif // OPM_WELL_FILTER_CAKE_HEADER_INCLUDED
