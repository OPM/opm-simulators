/*
  Copyright Equinor ASA 2021, 2022.

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/wells/SegmentState.hpp>

#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace {

std::vector<int> make_segment_number(const Opm::WellSegments& segments)
{
    std::vector<int> segment_number;
    segment_number.reserve(segments.size());

    std::transform(segments.begin(), segments.end(),
                   std::back_inserter(segment_number),
        [](const Opm::Segment& segment)
    {
        return segment.segmentNumber();
    });

    return segment_number;
}

} // Anonymous namespace

namespace Opm
{

SegmentState::SegmentState(int num_phases, const WellSegments& segments)
    : rates                    (segments.size() * num_phases)
    , dissolved_gas_rate       (segments.size())
    , vaporized_oil_rate       (segments.size())
    , phase_resv_rates         (segments.size() * num_phases)
    , phase_velocity           (segments.size() * num_phases)
    , phase_holdup             (segments.size() * num_phases)
    , phase_viscosity          (segments.size() * num_phases)
    , pressure                 (segments.size())
    , pressure_drop_friction   (segments.size())
    , pressure_drop_hydrostatic(segments.size())
    , pressure_drop_accel      (segments.size())
    , m_segment_number         (make_segment_number(segments))
{}

double SegmentState::pressure_drop(std::size_t index) const {
    return this->pressure_drop_friction[index] + this->pressure_drop_hydrostatic[index] + this->pressure_drop_accel[index];
}

bool SegmentState::empty() const {
    return this->rates.empty();
}

std::size_t SegmentState::size() const {
    return this->pressure.size();
}

void SegmentState::scale_pressure(const double bhp) {
    if (this->empty())
        throw std::logic_error("Tried to pressure scale empty SegmentState");

    const auto pressure_change = bhp - this->pressure[0];

    std::transform(this->pressure.begin(),
                   this->pressure.end(),
                   this->pressure.begin(),
                   [pressure_change] (const double& p) { return p + pressure_change;});
}

const std::vector<int>& SegmentState::segment_number() const {
    return this->m_segment_number;
}

} // namespace Opm
