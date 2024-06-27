/*
  Copyright 2023 Equinor ASA.

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

#include <opm/simulators/wells/ParallelPAvgCalculator.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/input/eclipse/EclipseState/Grid/GridDims.hpp>

#include <opm/input/eclipse/Schedule/Well/PAvgCalculator.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <functional>

template<class Scalar>
Opm::ParallelPAvgCalculator<Scalar>::
ParallelPAvgCalculator(const Parallel::Communication& comm,
                       const GridDims&                cellIndexMap,
                       const WellConnections&         connections)
    : PAvgCalculator<Scalar> { cellIndexMap, connections }
    , comm_          { comm }
{}

template<class Scalar>
void Opm::ParallelPAvgCalculator<Scalar>::
collectGlobalContributions()
{
    auto collect = [this](Accumulator& accumulator)
    {
        auto avg = accumulator.getRunningAverages();

        this->comm_.get().sum(avg.data(), avg.size());

        accumulator.assignRunningAverages(avg);
    };

    collect(this->accumCTF_);
    collect(this->accumPV_);
}

template class Opm::ParallelPAvgCalculator<double>;

#if FLOW_INSTANTIATE_FLOAT
template class Opm::ParallelPAvgCalculator<float>;
#endif
