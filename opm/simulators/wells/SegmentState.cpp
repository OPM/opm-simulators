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
    , phase_density            (segments.size() * (num_phases + 2)) // +2 for mixture with and without exponents
    , pressure                 (segments.size())
    , pressure_drop_friction   (segments.size())
    , pressure_drop_hydrostatic(segments.size())
    , pressure_drop_accel      (segments.size())
    , m_segment_number         (make_segment_number(segments))
{}

SegmentState SegmentState::serializationTestObject()
{
    SegmentState result;
    result.rates = {1.0, 2.0};
    result.dissolved_gas_rate = {3.0, 4.0, 5.0};
    result.vaporized_oil_rate = {6.0};
    result.phase_resv_rates = {7.0, 8.0};
    result.phase_velocity = {9.0};
    result.phase_holdup = {10.0, 11.0};
    result.phase_viscosity = {12.0, 12.5};
    result.phase_density = {13.0, 13.5, 13.6, 13.75};
    result.pressure = {14.0, 15.0};
    result.pressure_drop_friction = {16.0};
    result.pressure_drop_hydrostatic = {17.0, 18.0};
    result.pressure_drop_accel = {19.0};
    result.m_segment_number = {20, 21};

    return result;
}

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

bool SegmentState::operator==(const SegmentState& rhs) const
{
    return this->rates == rhs.rates &&
           this->dissolved_gas_rate == rhs.dissolved_gas_rate &&
           this->vaporized_oil_rate == rhs.vaporized_oil_rate &&
           this->phase_resv_rates == rhs.phase_resv_rates &&
           this->phase_velocity == rhs.phase_velocity &&
           this->phase_holdup == rhs.phase_holdup &&
           this->phase_viscosity == rhs.phase_viscosity &&
           this->phase_density == rhs.phase_density &&
           this->pressure == rhs.pressure &&
           this->pressure_drop_friction == rhs.pressure_drop_friction &&
           this->pressure_drop_hydrostatic == rhs.pressure_drop_hydrostatic &&
           this->pressure_drop_accel == rhs.pressure_drop_accel &&
           this->m_segment_number == rhs.m_segment_number;
}

} // namespace Opm
