/*
  Copyright 2026 Equinor ASA
  This file is part of OPM, distributed under the GNU General Public License,
  version 3 or (at your option) any later version.
*/
#ifndef OPM_GPU_BLACKOIL_NEWTON_VALIDATION_HPP
#define OPM_GPU_BLACKOIL_NEWTON_VALIDATION_HPP

#include <opm/models/blackoil/blackoilnewtonmethodparams.hpp>

namespace Opm::gpuistl {

// Bounded validation-only cases use separate buffers and leave resident state
// untouched. This runs in the same translation unit as the production kernel.
template<class CpuTypeTag, class Bridge, class Problem, class Solution>
void validateBlackoilNewtonBranches(Bridge& bridge, const Problem& problem,
                                   const Solution& solution,
                                   const BlackoilNewtonParams<typename Bridge::Scalar>& params);

// One-shot validation diagnostic for a solubility rounding difference. The
// two evaluations receive identical scalar inputs and neither edits state.
template<class CpuTypeTag, class Bridge>
void diagnoseBlackoilNewtonRsw(Bridge& bridge, unsigned region,
                              typename Bridge::Scalar temperature,
                              typename Bridge::Scalar pressure);

} // namespace Opm::gpuistl
#endif
