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
 * \copydoc Opm::RichardsNewtonMethod
 */
#ifndef EWOMS_RICHARDS_NEWTON_METHOD_HH
#define EWOMS_RICHARDS_NEWTON_METHOD_HH

#include "richardsproperties.hh"

#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/common/fvector.hh>

namespace Opm {

/*!
 * \ingroup RichardsModel
 *
 * \brief A Richards model specific Newton method.
 */
template <class TypeTag>
class RichardsNewtonMethod : public GET_PROP_TYPE(TypeTag, DiscNewtonMethod)
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscNewtonMethod) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Linearizer) Linearizer;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { pressureWIdx = Indices::pressureWIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { liquidPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };
    enum { gasPhaseIdx = GET_PROP_VALUE(TypeTag, GasPhaseIndex) };

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    RichardsNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {}

protected:
    friend NewtonMethod<TypeTag>;
    friend ParentType;

    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual OPM_UNUSED)
    {
        // normal Newton-Raphson update
        nextValue = currentValue;
        nextValue -= update;

        // do not clamp anything after 4 iterations
        if (this->numIterations_ > 4)
            return;

        const auto& problem = this->simulator_.problem();

        // calculate the old wetting phase saturation
        const MaterialLawParams& matParams =
            problem.materialLawParams(globalDofIdx, /*timeIdx=*/0);

        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;

        // set the temperature
        Scalar T = problem.temperature(globalDofIdx, /*timeIdx=*/0);
        fs.setTemperature(T);

        /////////
        // calculate the phase pressures of the previous iteration
        /////////

        // first, we have to find the minimum capillary pressure
        // (i.e. Sw = 0)
        fs.setSaturation(liquidPhaseIdx, 1.0);
        fs.setSaturation(gasPhaseIdx, 0.0);
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // non-wetting pressure can be larger than the
        // reference pressure if the medium is fully
        // saturated by the wetting phase
        Scalar pWOld = currentValue[pressureWIdx];
        Scalar pNOld =
            std::max(problem.referencePressure(globalDofIdx, /*timeIdx=*/0),
                     pWOld + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

        /////////
        // find the saturations of the previous iteration
        /////////
        fs.setPressure(liquidPhaseIdx, pWOld);
        fs.setPressure(gasPhaseIdx, pNOld);

        PhaseVector satOld;
        MaterialLaw::saturations(satOld, matParams, fs);
        satOld[liquidPhaseIdx] = std::max<Scalar>(0.0, satOld[liquidPhaseIdx]);

        /////////
        // find the wetting phase pressures which
        // corrospond to a 20% increase and a 20% decrease
        // of the wetting saturation
        /////////
        fs.setSaturation(liquidPhaseIdx, satOld[liquidPhaseIdx] - 0.2);
        fs.setSaturation(gasPhaseIdx, 1.0 - (satOld[liquidPhaseIdx] - 0.2));
        MaterialLaw::capillaryPressures(pC, matParams, fs);
        Scalar pwMin = pNOld - (pC[gasPhaseIdx] - pC[liquidPhaseIdx]);

        fs.setSaturation(liquidPhaseIdx, satOld[liquidPhaseIdx] + 0.2);
        fs.setSaturation(gasPhaseIdx, 1.0 - (satOld[liquidPhaseIdx] + 0.2));
        MaterialLaw::capillaryPressures(pC, matParams, fs);
        Scalar pwMax = pNOld - (pC[gasPhaseIdx] - pC[liquidPhaseIdx]);

        /////////
        // clamp the result to the minimum and the maximum
        // pressures we just calculated
        /////////
        Scalar pW = nextValue[pressureWIdx];
        pW = std::max(pwMin, std::min(pW, pwMax));
        nextValue[pressureWIdx] = pW;
    }
};
} // namespace Opm

#endif
