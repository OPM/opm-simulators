/*
  Copyright 2026 Equinor ASA
  This file is part of OPM, distributed under the GNU General Public License,
  version 3 or (at your option) any later version.
*/
#ifndef OPM_GPU_BLACKOIL_NEWTON_UPDATE_HPP
#define OPM_GPU_BLACKOIL_NEWTON_UPDATE_HPP

#include <opm/models/blackoil/blackoilnewtonupdate.hpp>

namespace Opm::gpuistl {

// Implemented in the Newton-only translation unit so its entire device call
// graph can preserve CPU arithmetic without changing the property kernel.
// The caller orders solver completion, reads compact status, then commits.
template<class DeviceTypeTag, class Bridge>
void launchBlackoilNewtonUpdate(Bridge& bridge,
                               const BlackoilNewtonParams<typename Bridge::Scalar>& params,
                               typename Bridge::Scalar relaxation,
                               bool useSOR, bool stabilize);

} // namespace Opm::gpuistl
#endif
