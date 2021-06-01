/*
  Copyright Equinor ASA 2021

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
#include <algorithm>
#include <stdexcept>

#include <opm/simulators/wells/SegmentState.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellConnections.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/WellSegments.hpp>

namespace Opm
{

SegmentState::SegmentState(int num_phases, const WellSegments& segments) :
    rates(segments.size() * num_phases),
    pressure(segments.size()),
    pressure_drop_friction(segments.size()),
    pressure_drop_hydrostatic(segments.size()),
    pressure_drop_accel(segments.size())
{

}

double SegmentState::pressure_drop(std::size_t index) const {
    return this->pressure_drop_friction[index] + this->pressure_drop_hydrostatic[index] + this->pressure_drop_accel[index];
}


bool SegmentState::empty() const {
    return this->rates.empty();
}

void SegmentState::scale_pressure(double bhp) {
    if (this->empty())
        throw std::logic_error("Tried to pressure scale empty SegmentState");

    auto scale_factor = bhp / this->pressure[0];
    std::transform(this->pressure.begin(),
                   this->pressure.end(),
                   this->pressure.begin(),
                   [scale_factor] (const double& p) { return p*scale_factor;});
}

}
