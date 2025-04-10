// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/simulators/flow/TracerContainer.hpp>

#include <opm/input/eclipse/EclipseState/TracerConfig.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <algorithm>
#include <utility>

namespace Opm {

template<class FluidSystem>
void TracerContainer<FluidSystem>::
allocate(const unsigned bufferSize,
         const TracerConfig& tracers)
{
    if (!tracers.empty()) {
        allocated_ = true;
        freeConcentrations_.resize(tracers.size());
        solConcentrations_.resize(tracers.size());
        std::for_each(tracers.begin(), tracers.end(),
                      [idx = 0, bufferSize, this](const auto& tracer) mutable
                      {
                          freeConcentrations_[idx].resize(bufferSize, 0.0);
                          if (((tracer.phase == Phase::GAS && FluidSystem::enableDissolvedGas()) ||
                               (tracer.phase == Phase::OIL && FluidSystem::enableVaporizedOil())) &&
                              (tracer.solution_concentration.has_value() ||
                               tracer.solution_tvdp.has_value()))
                          {
                              solConcentrations_[idx].resize(bufferSize, 0.0);
                          }
                          ++idx;
                      });
    }
}

template<class FluidSystem>
void TracerContainer<FluidSystem>::
assignFreeConcentrations(const unsigned globalDofIdx,
                         const AssignFunction& concentration)
{
    std::for_each(freeConcentrations_.begin(), freeConcentrations_.end(),
                  [globalDofIdx, idx = 0, &concentration](auto& tracer) mutable
                  {
                      if (!tracer.empty()) {
                          tracer[globalDofIdx] = concentration(idx);
                      }
                      ++idx;
                  });
}

template<class FluidSystem>
void TracerContainer<FluidSystem>::
assignSolConcentrations(const unsigned globalDofIdx,
                         const AssignFunction& concentration)
{
    std::for_each(solConcentrations_.begin(), solConcentrations_.end(),
                  [globalDofIdx, idx = 0, &concentration](auto& tracer) mutable
                  {
                      if (!tracer.empty()) {
                          tracer[globalDofIdx] = concentration(idx);
                      }
                      ++idx;
                  });
}

template<class FluidSystem>
void TracerContainer<FluidSystem>::
outputRestart(data::Solution& sol,
              const TracerConfig& tracers)
{
    if (!this->allocated_) {
        return;
    }

    std::for_each(tracers.begin(), tracers.end(),
                    [idx = 0, &sol, this](const auto& tracer) mutable
                    {
                        sol.insert(tracer.fname(),
                                   UnitSystem::measure::identity,
                                   std::move(freeConcentrations_[idx]),
                                   data::TargetType::RESTART_TRACER_SOLUTION);

                        if (!solConcentrations_[idx].empty()) {
                            sol.insert(tracer.sname(),
                                       UnitSystem::measure::identity,
                                       std::move(solConcentrations_[idx]),
                                       data::TargetType::RESTART_TRACER_SOLUTION);
                        }
                        ++idx;
                    });

    this->allocated_ = false;
}

template<class T> using FS = BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>;

#define INSTANTIATE_TYPE(T) \
    template class TracerContainer<FS<T>>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

#define INSTANTIATE_COMP_THREEPHASE(NUM) \
    template<class T> using FS##NUM = GenericOilGasWaterFluidSystem<T, NUM, true>; \
    template class TracerContainer<FS##NUM<double>>;

#define INSTANTIATE_COMP_TWOPHASE(NUM) \
    template<class T> using GFS##NUM = GenericOilGasWaterFluidSystem<T, NUM, false>; \
    template class TracerContainer<GFS##NUM<double>>;

#define INSTANTIATE_COMP(NUM) \
    INSTANTIATE_COMP_THREEPHASE(NUM) \
    INSTANTIATE_COMP_TWOPHASE(NUM)

INSTANTIATE_COMP_THREEPHASE(0)  // \Note: to register the parameter ForceDisableFluidInPlaceOutput
INSTANTIATE_COMP(2)
INSTANTIATE_COMP(3)
INSTANTIATE_COMP(4)
INSTANTIATE_COMP(5)
INSTANTIATE_COMP(6)
INSTANTIATE_COMP(7)

} // namespace Opm
