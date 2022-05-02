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
#ifndef EWOMS_FLASH_NEWTON_METHOD_HH
#define EWOMS_FLASH_NEWTON_METHOD_HH

#include "flashproperties.hh"

#include <opm/models/nonlinear/newtonmethod.hh>

#include <opm/common/Exceptions.hpp>

#include <algorithm>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup FlashModel
 *
 * \brief A Newton solver specific to the NCP model.
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

public:
    /*!
     * \copydoc FvBaseNewtonMethod::FvBaseNewtonMethod(Problem& )
     */
    FlashNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {}

protected:
    friend ParentType;
    friend NewtonMethod<TypeTag>;

    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        // normal Newton-Raphson update
        nextValue = currentValue;
        nextValue -= update;

        ////
        // Pressure updates
        ////
        // limit pressure reference change to 20% of the total value per iteration
        clampValue_(nextValue[pressure0Idx],
                    currentValue[pressure0Idx]*0.8,
                    currentValue[pressure0Idx]*1.2);

        ////
        // z updates
        ////
        // restrict update to at most 0.1
        Scalar maxDeltaZ = 0.0;  // in update vector
        Scalar sumDeltaZ = 0.0; // changes in last component (not in update vector)
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
            maxDeltaZ = std::max(std::abs(update[z0Idx + compIdx]), maxDeltaZ);
            sumDeltaZ += update[z0Idx + compIdx];
        }
        maxDeltaZ = std::max(std::abs(-sumDeltaZ), maxDeltaZ);

        // if max. update is above 0.1, restrict that one to 0.1 and adjust the rest accordingly (s.t. last comp. update is sum of the changes in update vector)
        if (maxDeltaZ > 0.1) {
            Scalar alpha = 0.2/maxDeltaZ;
            for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
                nextValue[z0Idx + compIdx] = currentValue[z0Idx + compIdx] - alpha*update[z0Idx + compIdx];
            }
        }

        // ensure that z-values are less than tol or more than 1-tol
        Scalar tol = 1e-8;
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
            clampValue_(nextValue[z0Idx + compIdx], tol, 1-tol);
        }

        // normalize s.t. z sums to 1.0
        // Scalar lastZ = 1.0;
        // Scalar sumZ = 0.0;
        // for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
        //     lastZ -= nextValue[z0Idx + compIdx];
        //     sumZ += nextValue[z0Idx + compIdx];
        // }
        // sumZ += lastZ;
        // for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
        //     nextValue[z0Idx + compIdx] /= sumZ;
        // }
    }
private:
    void clampValue_(Scalar& val, Scalar minVal, Scalar maxVal) const
    { val = std::max(minVal, std::min(val, maxVal)); }

};  // class FlashNewtonMethod
} // namespace Opm
#endif