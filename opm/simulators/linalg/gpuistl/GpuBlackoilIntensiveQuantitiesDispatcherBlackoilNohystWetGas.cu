/*
  Copyright 2026 Equinor ASA

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
/*!
 * \file
 *
 * GPU intensive-quantities dispatcher configuration for three-phase black-oil
 * decks with live oil and wet gas (e.g. Norne), run by \c flow with
 * \c FlowProblemNohyst.
 */
#include <config.h>

// Tricks for older versions of ROCm that do not support alugrid or dune-fem
// (mirrors test_blackoilintensivequantities_gpu.cu).
#ifdef HAVE_DUNE_ALUGRID
#undef HAVE_DUNE_ALUGRID
#endif
#define HAVE_DUNE_ALUGRID 0
#ifdef HAVE_DUNE_FEM
#undef HAVE_DUNE_FEM
#endif
#define HAVE_DUNE_FEM 0

#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher_impl.hpp>

#include <opm/material/fluidsystems/blackoilpvt/BlackOilPvtApproaches.hpp>

#include <opm/simulators/flow/FlowProblemNohystTypeTag.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowBlackoilNohystTypeTag.hpp>

namespace Opm::gpuistl {

template <>
struct GpuDispatcherTraits<Properties::TTag::FlowProblemNohyst, LiveOilWetGasPvtApproaches>
    : dispatcher_detail::BlackoilNohystDispatcherTraits<
          Properties::TTag::FlowBlackoilNohystWetGasGPU>
{
    static constexpr const char* name = "three-phase black-oil: live oil and wet gas";
};

// Used by flow for black-oil decks without HYSTER.
template std::unique_ptr<
    GpuBlackoilIntensiveQuantitiesDispatcher<Properties::TTag::FlowProblemNohyst>::Impl>
GpuBlackoilIntensiveQuantitiesDispatcher<Properties::TTag::FlowProblemNohyst>::
    makeImpl<LiveOilWetGasPvtApproaches>(const Problem&);

} // namespace Opm::gpuistl
