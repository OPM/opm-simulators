/*
  Copyright 2022-2023 SINTEF AS

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

#ifndef OPM_CUISTL_PRECONDIDTIONER_SHOULD_CALL_POST_PRE_HPP
#define OPM_CUISTL_PRECONDIDTIONER_SHOULD_CALL_POST_PRE_HPP

#include <opm/simulators/linalg/cuistl/detail/has_function.hpp>

namespace Opm::cuistl::detail
{

//! @brief Tests (compile time) if the preconditioner type needs to call pre() before a call to apply()
//!
//! @note This is mostly used to avoid unneeded copying back and front to the GPU, as well
//! as avoiding communication.
template <class PreconditionerType>
constexpr bool
shouldCallPreconditionerPre()
{
    if constexpr (has_should_call_pre<PreconditionerType>::value) {
        return PreconditionerType::shouldCallPre();
    } else {
        // If the preconditioner type does not have a way of signalling the need
        // to call pre(), we should always call it. This is because it could be an
        // old design.
        return true;
    }
}

//! @brief Tests (compile time) if the preconditioner type needs to call post() after a call to apply(...)
//!
//! @note This is mostly used to avoid unneeded copying back and front to the GPU, as well
//! as avoiding communication.
template <class PreconditionerType>
constexpr bool
shouldCallPreconditionerPost()
{
    if constexpr (has_should_call_post<PreconditionerType>::value) {
        return PreconditionerType::shouldCallPost();
    } else {
        // If the preconditioner type does not have a way of signalling the need
        // to call post(), we should always call it. This is because it could be an
        // old design.
        return true;
    }
}
} // namespace Opm::cuistl::detail
#endif
