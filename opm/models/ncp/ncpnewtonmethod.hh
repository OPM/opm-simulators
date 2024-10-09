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
 * \copydoc Opm::NcpNewtonMethod
 */
#ifndef EWOMS_NCP_NEWTON_METHOD_HH
#define EWOMS_NCP_NEWTON_METHOD_HH

#include "ncpproperties.hh"

#include <opm/common/Exceptions.hpp>

#include <opm/models/nonlinear/newtonmethod.hh>

#include <algorithm>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup NcpModel
 *
 * \brief A Newton solver specific to the NCP model.
 */
template <class TypeTag>
class NcpNewtonMethod : public GetPropType<TypeTag, Properties::DiscNewtonMethod>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscNewtonMethod>;

    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { fugacity0Idx = Indices::fugacity0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { ncp0EqIdx = Indices::ncp0EqIdx };

public:
    /*!
     * \copydoc FvBaseNewtonMethod::FvBaseNewtonMethod(Problem& )
     */
    NcpNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {}

protected:
    friend ParentType;
    friend NewtonMethod<TypeTag>;

    void preSolve_(const SolutionVector&,
                   const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = this->model().linearizer().constraintsMap();
        this->lastError_ = this->error_;

        // calculate the error as the maximum weighted tolerance of
        // the solution's residual
        this->error_ = 0;
        for (unsigned dofIdx = 0; dofIdx < currentResidual.size(); ++dofIdx) {
            // do not consider auxiliary DOFs for the error
            if (dofIdx >= this->model().numGridDof() || this->model().dofTotalVolume(dofIdx) <= 0.0)
                continue;

            // also do not consider DOFs which are constraint
            if (this->enableConstraints_()) {
                if (constraintsMap.count(dofIdx) > 0)
                    continue;
            }

            const auto& r = currentResidual[dofIdx];
            for (unsigned eqIdx = 0; eqIdx < r.size(); ++eqIdx) {
                if (ncp0EqIdx <= eqIdx && eqIdx < Indices::ncp0EqIdx + numPhases)
                    continue;
                this->error_ =
                    std::max(std::abs(r[eqIdx]*this->model().eqWeight(dofIdx, eqIdx)),
                             this->error_);
            }
        }

        // take the other processes into account
        this->error_ = this->comm_.max(this->error_);

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (this->error_ > Parameters::Get<Parameters::NewtonMaxError<Scalar>>())
            throw Opm::NumericalProblem("Newton: Error "+std::to_string(double(this->error_))+
                                        + " is larger than maximum allowed error of "
                                        + std::to_string(Parameters::Get<Parameters::NewtonMaxError<Scalar>>()));
    }

    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector&)
    {
        // normal Newton-Raphson update
        nextValue = currentValue;
        nextValue -= update;

        ////
        // put crash barriers along the update path
        ////

        // saturations: limit the change of any saturation to at most 20%
        Scalar sumSatDelta = 0.0;
        Scalar maxSatDelta = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            maxSatDelta = std::max(std::abs(update[saturation0Idx + phaseIdx]),
                                   maxSatDelta);
            sumSatDelta += update[saturation0Idx + phaseIdx];
        }
        maxSatDelta = std::max(std::abs(- sumSatDelta), maxSatDelta);

        if (maxSatDelta > 0.2) {
            Scalar alpha = 0.2/maxSatDelta;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
                nextValue[saturation0Idx + phaseIdx] =
                    currentValue[saturation0Idx + phaseIdx]
                    - alpha*update[saturation0Idx + phaseIdx];
            }
        }

        // limit pressure reference change to 20% of the total value per iteration
        clampValue_(nextValue[pressure0Idx],
                    currentValue[pressure0Idx]*0.8,
                    currentValue[pressure0Idx]*1.2);

        // fugacities
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar& val = nextValue[fugacity0Idx + compIdx];
            Scalar oldVal = currentValue[fugacity0Idx + compIdx];

            // get the minimum activity coefficient for the component (i.e., the activity
            // coefficient of the phase for which the component has the highest affinity)
            Scalar minPhi = this->problem().model().minActivityCoeff(globalDofIdx, compIdx);
            // Make sure that the activity coefficient does not get too small.
            minPhi = std::max(0.001*currentValue[pressure0Idx], minPhi);

            // allow the mole fraction of the component to change at most 70% in any
            // phase (assuming composition independent fugacity coefficients).
            Scalar maxDelta = 0.7 * minPhi;
            clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);

            // make sure that fugacities do not become negative
            val = std::max(val, 0.0);
        }

        // do not become grossly unphysical in a single iteration for the first few
        // iterations of a time step
        if (this->numIterations_ < 3) {
            // fugacities
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar& val = nextValue[fugacity0Idx + compIdx];
                Scalar oldVal = currentValue[fugacity0Idx + compIdx];
                Scalar minPhi = this->problem().model().minActivityCoeff(globalDofIdx, compIdx);
                if (oldVal < 1.0*minPhi && val > 1.0*minPhi)
                    val = 1.0*minPhi;
                else if (oldVal > 0.0 && val < 0.0)
                    val = 0.0;
            }

            // saturations
            for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
                Scalar& val = nextValue[saturation0Idx + phaseIdx];
                Scalar oldVal = currentValue[saturation0Idx + phaseIdx];
                if (oldVal < 1.0 && val > 1.0)
                    val = 1.0;
                else if (oldVal > 0.0 && val < 0.0)
                    val = 0.0;
            }
        }
    }

private:
    void clampValue_(Scalar& val, Scalar minVal, Scalar maxVal) const
    { val = std::max(minVal, std::min(val, maxVal)); }
};
} // namespace Opm

#endif
