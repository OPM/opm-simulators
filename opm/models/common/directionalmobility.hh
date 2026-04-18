// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 Equinor ASA.

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
 * \brief This file contains definitions related to directional mobilities
 */
#ifndef OPM_MODELS_DIRECTIONAL_MOBILITY_HH
#define OPM_MODELS_DIRECTIONAL_MOBILITY_HH

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <array>
#include <stdexcept>

namespace Opm {

template <class TypeTag>
class DirectionalMobility
{
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    using array_type = std::array<Evaluation,numPhases>;

    DirectionalMobility() = default;

    DirectionalMobility(const array_type& mX,
                        const array_type& mY,
                        const array_type& mZ)
        : mobility_{mX, mY, mZ}
    {}

    OPM_HOST_DEVICE const array_type& getArray(unsigned index) const
    {
        if (index > 2) {
            OPM_THROW(std::runtime_error, "Unexpected mobility array index");
        }
        return mobility_[index];
    }

    OPM_HOST_DEVICE array_type& getArray(unsigned index)
    {
        return const_cast<array_type&>(std::as_const(*this).getArray(index));
    }

private:
    std::array<array_type,3> mobility_{};
};

} // namespace Opm

#endif
