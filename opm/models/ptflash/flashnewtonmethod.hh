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

#include <algorithm>
#include <cmath>

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
        // normal Newton-Raphson update
        nextValue = currentValue;
        nextValue -= update;

        ////
        // Pressure updates
        ////
        // limit pressure reference change relative to the total value per iteration
        constexpr Scalar max_percent_change = 0.2;
        constexpr Scalar upper_bound = 1. + max_percent_change;
        constexpr Scalar lower_bound = 1. - max_percent_change;
        nextValue[pressure0Idx] = std::clamp(nextValue[pressure0Idx],
                                             currentValue[pressure0Idx] * lower_bound,
                                             currentValue[pressure0Idx] * upper_bound);

        ////
        // z updates
        ////
        // restrict update
        Scalar maxDeltaZ = 0.0;  // in update vector
        Scalar sumDeltaZ = 0.0; // changes in last component (not in update vector)
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
            maxDeltaZ = std::max(std::abs(update[z0Idx + compIdx]), maxDeltaZ);
            sumDeltaZ += update[z0Idx + compIdx];
        }
        maxDeltaZ = std::max(std::abs(sumDeltaZ), maxDeltaZ);

        // if max. update is above limit, restrict that one to limit and adjust the rest
        // accordingly (s.t. last comp. update is sum of the changes in update vector)
        // \Note: original code uses 0.1, while 0.1 looks like having problem make it converged.
        // So there is some more to investigate here
        constexpr Scalar deltaz_limit = 0.2;
        if (maxDeltaZ > deltaz_limit) {
            const Scalar alpha = deltaz_limit / maxDeltaZ;
            for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
                nextValue[z0Idx + compIdx] = currentValue[z0Idx + compIdx] - alpha * update[z0Idx + compIdx];
            }
        }

        // ensure that z-values are less than tol or more than 1-tol
        constexpr Scalar tol = 1e-8;
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
           nextValue[z0Idx + compIdx] = std::clamp(nextValue[z0Idx + compIdx], tol, 1-tol);
        }

        if constexpr (waterEnabled) {
            // limit change in water saturation
            constexpr Scalar dSwMax = 0.2;
            if (update[Indices::water0Idx] > dSwMax) {
                nextValue[Indices::water0Idx] = currentValue[Indices::water0Idx] - dSwMax;
            }
        }
    }
};  // class FlashNewtonMethod

} // namespace Opm

#endif
