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
 * \copydoc Opm::FlashNewtonMethod
 */
#ifndef OPM_FLASH_NEWTON_METHOD_HH
#define OPM_FLASH_NEWTON_METHOD_HH

#include <opm/common/Exceptions.hpp>

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/nonlinear/newtonmethod.hh>
#include <opm/models/ptflash/flashcompositionstep.hh>

#include <algorithm>
#include <array>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup FlashModel
 *
 * \brief A Newton solver specific to the PTFlash model.
 */
template <class TypeTag>
class FlashNewtonMethod : public GetPropType<TypeTag, Properties::DiscNewtonMethod>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscNewtonMethod>;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    enum { pressure0Idx = Indices::pressure0Idx };
    enum { z0Idx = Indices::z0Idx };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    static constexpr bool waterEnabled = Indices::waterEnabled;

public:
    /*!
     * \copydoc FvBaseNewtonMethod::FvBaseNewtonMethod(Problem& )
     */
    explicit FlashNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {}

protected:
    friend ParentType;
    friend NewtonMethod<TypeTag>;

    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned /* globalDofIdx */,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& /* currentResidual */)
    {
        // nextValue may alias currentValue. Preserve the pre-update value because the
        // limiters below must chop relative to it after nextValue has been updated.
        const PrimaryVariables priVarsOld = currentValue;

        // normal Newton-Raphson update
        nextValue = priVarsOld;
        nextValue -= update;

        ////
        // Pressure updates
        ////
        // limit pressure reference change relative to the total value per iteration
        constexpr Scalar max_percent_change = 0.2;
        constexpr Scalar upper_bound = 1. + max_percent_change;
        constexpr Scalar lower_bound = 1. - max_percent_change;
        nextValue[pressure0Idx] = std::clamp(nextValue[pressure0Idx],
                                             priVarsOld[pressure0Idx] * lower_bound,
                                             priVarsOld[pressure0Idx] * upper_bound);

        ////
        // Composition and water saturation updates
        ////
        // Steps are new minus old values. The last fraction is dependent: it and its step
        // complete the sums of z and dz to one and zero.
        std::array<Scalar, numComponents> z{};
        std::array<Scalar, numComponents> dz{};
        z.back() = 1.0;
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
            z[compIdx] = priVarsOld[z0Idx + compIdx];
            dz[compIdx] = -update[z0Idx + compIdx];
            z.back() -= z[compIdx];
            dz.back() -= dz[compIdx];
        }

        Scalar sw = 0.0;
        Scalar dSw = 0.0;
        if constexpr (waterEnabled) {
            sw = priVarsOld[Indices::water0Idx];
            dSw = -update[Indices::water0Idx];
        }

        constexpr Scalar maxAmountChange = 0.2;
        applyFlashCompositionStep(z, sw, dz, dSw,
                                  IntensiveQuantities::compositionFloor,
                                  IntensiveQuantities::hydrocarbonFloor,
                                  maxAmountChange);

        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
            nextValue[z0Idx + compIdx] = z[compIdx];
        }
        if constexpr (waterEnabled) {
            nextValue[Indices::water0Idx] = sw;
        }
    }
};  // class FlashNewtonMethod

} // namespace Opm

#endif
