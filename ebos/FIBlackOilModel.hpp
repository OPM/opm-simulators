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
/*!
 * \file
 *
 * \copydoc Opm::FIBlackOilModel
 */
#ifndef FI_BLACK_OIL_MODEL_HPP
#define FI_BLACK_OIL_MODEL_HPP

#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/common/ErrorMacros.hpp>

#include <stdexcept>

namespace Opm{
    template<typename TypeTag>
    class FIBlackOilModel: public BlackOilModel<TypeTag>{
        using Parent = BlackOilModel<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    public:
        FIBlackOilModel(Simulator& simulator): BlackOilModel<TypeTag>(simulator){
        }

        // standard flow
        const IntensiveQuantities& intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const{
            if (!this->enableIntensiveQuantityCache_){
                OPM_THROW(std::logic_error, "Run without intensive quantites not enabled: Use --enable-intensive-quantity=true");
            }
            const auto* intquant = this->cachedIntensiveQuantities(globalIdx, timeIdx);
            if(!intquant){
                OPM_THROW(std::logic_error, "Intensive quantites need to be updated in code");
            }
            return *intquant;
        }

    };
} // namespace Opm
#endif // FI_BLACK_OIL_MODEL_HPP
