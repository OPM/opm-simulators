/*
  Copyright 2026 Equinor ASA.

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

#ifndef OPM_WELL_RESERVOIR_RESPONSE_HEADER_INCLUDED
#define OPM_WELL_RESERVOIR_RESPONSE_HEADER_INCLUDED

#include <optional>
#include <vector>

namespace Opm
{

/// Response of the pressures in the cells of a well's connections to the flow
/// of the well, within a time step.
///
/// The flow of a well lowers the pressures in its cells within the same time
/// step, so a well that can flow at the cell pressures of the stopped well may
/// not be able to flow at those its own flow leaves.  The response is measured
/// between two consecutive Newton iterations in which the well is stopped in
/// one and flowing in the other: the change of the cell pressures divided by
/// the drawdown of the flowing well.  With the inflow linear in the bhp at
/// fixed cell pressures,
///
///     q_c = a_c - b_c * bhp
///
/// for each component c, pressures and drawdowns are averages over the
/// components weighted by b_c.
template <class Scalar>
class WellReservoirResponse
{
public:
    /// Record the well at the start of a Newton iteration.
    ///
    /// \param[in] iteration Newton iteration in the time step.
    /// \param[in] stopped Whether the well is stopped.
    /// \param[in] ipr_a a_c at the current cell pressures.
    /// \param[in] ipr_b b_c at the current cell pressures.
    /// \param[in] rates Surface rate per component, production positive.
    void update(int iteration,
                bool stopped,
                const std::vector<Scalar>& ipr_a,
                const std::vector<Scalar>& ipr_b,
                const std::vector<Scalar>& rates);

    /// Change of the cell pressures per unit drawdown of the well, once measured
    /// and while the well has an inflow.
    std::optional<Scalar> ratio() const
    {
        return ratio_;
    }

    /// Drop of the cell pressures from those of the last update when the well
    /// flows at a bhp, negative where they rise.  Lowering all cell pressures
    /// by this drop changes the inflow as raising the bhp by it does.
    /// Requires ratio().
    Scalar pressureDrop(Scalar bhp) const;

    /// Rise from the cell pressures of the last update to those of the stopped
    /// well.  Requires ratio().
    Scalar stoppedPressureRise() const;

private:
    struct Sample {
        int iteration;
        bool stopped;
        Scalar pressure;
        Scalar drawdown;
    };

    std::optional<Sample> last_ {};
    std::optional<Scalar> ratio_ {};
};

} // namespace Opm

#endif // OPM_WELL_RESERVOIR_RESPONSE_HEADER_INCLUDED
