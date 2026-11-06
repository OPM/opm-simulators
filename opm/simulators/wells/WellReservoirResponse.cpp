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

#include <config.h>

#include <opm/simulators/wells/WellReservoirResponse.hpp>

#include <algorithm>
#include <cassert>
#include <numeric>

namespace Opm
{

template <class Scalar>
void
WellReservoirResponse<Scalar>::update(const int iteration,
                                      const bool stopped,
                                      const std::vector<Scalar>& ipr_a,
                                      const std::vector<Scalar>& ipr_b,
                                      const std::vector<Scalar>& rates)
{
    const Scalar sum_b = std::accumulate(ipr_b.begin(), ipr_b.end(), Scalar {0});
    if (sum_b <= 0) {
        // Without inflow there are no cell pressures to predict from.
        last_.reset();
        ratio_.reset();
        return;
    }
    const Scalar pressure = std::accumulate(ipr_a.begin(), ipr_a.end(), Scalar {0}) / sum_b;
    const Scalar drawdown = stopped
        ? Scalar {0}
        : std::max(Scalar {0}, std::accumulate(rates.begin(), rates.end(), Scalar {0}) / sum_b);

    // At iteration 0 the cell pressures are still those of the previous time step. Later
    // changes of status are damped by the relaxation of oscillating Newton updates.
    if (!ratio_ && last_ && last_->iteration > 0 && last_->iteration + 1 == iteration
        && last_->stopped != stopped) {
        const Scalar rise = stopped ? pressure - last_->pressure : last_->pressure - pressure;
        const Scalar flowing_drawdown = stopped ? last_->drawdown : drawdown;
        if (rise > 0 && flowing_drawdown > 0) {
            ratio_ = rise / flowing_drawdown;
        }
    }
    last_ = Sample {iteration, stopped, pressure, drawdown};
}

template <class Scalar>
Scalar
WellReservoirResponse<Scalar>::pressureDrop(const Scalar bhp) const
{
    assert(ratio_ && last_);
    // With the drawdown d from the cell pressure p of the stopped well, the cells are at
    // p - ratio * d, which gives d = (p - bhp) / (1 + ratio).
    const Scalar stopped_pressure = last_->pressure + stoppedPressureRise();
    return *ratio_ / (1 + *ratio_) * (stopped_pressure - bhp) - stoppedPressureRise();
}

template <class Scalar>
Scalar
WellReservoirResponse<Scalar>::stoppedPressureRise() const
{
    assert(ratio_ && last_);
    return *ratio_ * last_->drawdown;
}

template class WellReservoirResponse<double>;

#if FLOW_INSTANTIATE_FLOAT
template class WellReservoirResponse<float>;
#endif

} // namespace Opm
