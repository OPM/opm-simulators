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

#ifndef OPM_SEGMENTSTATE_HEADER_INCLUDED
#define OPM_SEGMENTSTATE_HEADER_INCLUDED

#include <cstddef>
#include <vector>

namespace Opm
{
class WellSegments;
} // namespace Opm

namespace Opm
{

class SegmentState
{
public:
    SegmentState() = default;
    SegmentState(int num_phases, const WellSegments& segments);

    static SegmentState serializationTestObject();

    double pressure_drop(std::size_t index) const;
    bool empty() const;
    void scale_pressure(double bhp);

    const std::vector<int>& segment_number() const;
    std::size_t size() const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(rates);
        serializer(dissolved_gas_rate);
        serializer(vaporized_oil_rate);
        serializer(phase_resv_rates);
        serializer(phase_velocity);
        serializer(phase_holdup);
        serializer(phase_viscosity);
        serializer(phase_density);
        serializer(pressure);
        serializer(pressure_drop_friction);
        serializer(pressure_drop_hydrostatic);
        serializer(pressure_drop_accel);
        serializer(m_segment_number);
    }

    bool operator==(const SegmentState&) const;

    std::vector<double> rates;
    std::vector<double> dissolved_gas_rate;
    std::vector<double> vaporized_oil_rate;

    /// Segment condition volume flow rates through segment (per phase)
    std::vector<double> phase_resv_rates;

    /// Segment condition flow velocity through segment (per phase)
    std::vector<double> phase_velocity;

    /// Segment condition holdup fractions through segment (per phase)
    std::vector<double> phase_holdup;

    /// Segment condition phase viscosities.
    std::vector<double> phase_viscosity;

    /// Segment condition phase densities.
    ///
    /// \code num_phases + 2 \endcode elements per segment, with the two
    /// additional elements being fluid mixture densities.  Fluid mixture
    /// densities are calculated with and without ICD-related flowing phase
    /// fraction exponents, respectively.  Mixture density *without* flowing
    /// phase fraction exponents is stored at offset \c num_phases, and
    /// mixture density *with* flowing phase fraction exponents is stored at
    /// offset \code num_phases + 1 \endcode.
    ///
    /// In other words, the phase and mixture densities are stored as
    ///
    ///   [{ p0, p1, ..., (np - 1), mixture, mixture_with_exponents },
    ///    { p0, p1, ..., (np - 1), mixture, mixture_with_exponents },
    ///    ...
    ///    { p0, p1, ..., (np - 1), mixture, mixture_with_exponents }]
    std::vector<double> phase_density;

    std::vector<double> pressure;
    std::vector<double> pressure_drop_friction;
    std::vector<double> pressure_drop_hydrostatic;
    std::vector<double> pressure_drop_accel;

private:
    std::vector<int>    m_segment_number;
};

} // namepace Opm

#endif  // OPM_SEGMENTSTATE_HEADER_INCLUDED
