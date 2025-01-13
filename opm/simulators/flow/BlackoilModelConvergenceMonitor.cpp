/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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
#include <opm/simulators/flow/BlackoilModelConvergenceMonitor.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace Opm {

template<class Scalar>
BlackoilModelConvergenceMonitor<Scalar>::
BlackoilModelConvergenceMonitor(const MonitorParams& param)
    : param_(param)
    , prev_distance_(std::numeric_limits<double>::infinity())
    , prev_above_tolerance_(0)
{
}

template<class Scalar>
void BlackoilModelConvergenceMonitor<Scalar>::
checkPenaltyCard(ConvergenceReport& report, int iteration)
{
    const auto& current_metrics = report.reservoirConvergence();
    auto distances = std::vector<double>(current_metrics.size(), 0.0);
    int current_above_tolerance = 0;

    for (size_t i = 0; i < current_metrics.size(); ++i) {
        distances[i] = std::max(std::log10(current_metrics[i].value() /
                                           current_metrics[i].tolerance()), 0.0);
            // Count number of metrics above tolerance
            if (current_metrics[i].value() > current_metrics[i].tolerance()) {
                current_above_tolerance++;
            }
        }

    // use L1 norm of the distances vector
    double current_distance = std::accumulate(distances.begin(), distances.end(), 0.0);

    if (iteration > 0) {
        // Add penalty if number of metrics above tolerance has increased
        if (current_above_tolerance > prev_above_tolerance_) {
            report.addNonConvergedPenalty();
        }

        if (current_distance > param_.decay_factor_ * prev_distance_) {
            report.addDistanceDecayPenalty();
        }
    }

    prev_distance_ = current_distance;
    prev_above_tolerance_ = current_above_tolerance;

    if (report.wellFailures().size() > 0) {
        report.addLargeWellResidualsPenalty();
    }

    total_penaltyCard_ += report.getPenaltyCard();

    if (param_.enabled_ && (total_penaltyCard_.total() > param_.cutoff_)) {
        report.setReservoirFailed(
            {ConvergenceReport::ReservoirFailure::Type::ConvergenceMonitorFailure,
             ConvergenceReport::Severity::ConvergenceMonitorFailure,
             -1}); // -1 indicates it's not specific to any component
    }
}

template<class Scalar>
void BlackoilModelConvergenceMonitor<Scalar>::
reset()
{
    total_penaltyCard_.reset();
    prev_above_tolerance_ = 0;
    prev_distance_ = std::numeric_limits<double>::infinity();
}

template class BlackoilModelConvergenceMonitor<double>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilModelConvergenceMonitor<float>;
#endif

} // namespace Opm
