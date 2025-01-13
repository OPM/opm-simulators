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
#include <opm/models/utils/propertysystem.hh>
#include <opm/material/densead/Evaluation.hpp>

#include <array>
#include <stdexcept>

namespace Opm {
template <class TypeTag, class Evaluation>
struct DirectionalMobility {
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    // TODO: This (line below) did not work. I get error: ‘Evaluation’ is not a member of ‘Opm::Properties’
    //  when compiling the tracer model (eclgenerictracermodel.cc). 
    //  QuickFix: I am adding Evaluation as a class template parameter..
    //using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

    using array_type = std::array<Evaluation,numPhases>;
    DirectionalMobility(const DirectionalMobility& other)
        : mobilityX_{other.mobilityX_}, mobilityY_{other.mobilityY_}, mobilityZ_{other.mobilityZ_} {}
    DirectionalMobility(const array_type& mX, const array_type& mY, const array_type& mZ)
        : mobilityX_{mX}, mobilityY_{mY}, mobilityZ_{mZ} {}
    DirectionalMobility() : mobilityX_{}, mobilityY_{}, mobilityZ_{} {}
    array_type& getArray(int index) {
        switch(index) {
            case 0:
                return mobilityX_;
            case 1:
                return mobilityY_;
            case 2:
                return mobilityZ_;
            default:
                throw std::runtime_error("Unexpected mobility array index");
        }
    }
    array_type mobilityX_;
    array_type mobilityY_;
    array_type mobilityZ_;
};
} // namespace Opm
#endif
