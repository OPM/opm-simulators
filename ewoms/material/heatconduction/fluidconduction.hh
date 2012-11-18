// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Ewoms::FluidHeatConduction
 */
#ifndef EWOMS_FLUID_HEAT_CONDUCTION_HH
#define EWOMS_FLUID_HEAT_CONDUCTION_HH

#include "fluidconductionparams.hh"

#include <ewoms/common/spline.hh>
#include <algorithm>

namespace Ewoms
{
/*!
 * \ingroup material
 *
 * \brief Implements a heat conduction law which just takes the conductivity of a given fluid phase.
 */
template <class FluidSystem,
          class ScalarT,
          int phaseIdx,
          class ParamsT = FluidHeatConductionParams<ScalarT> >
class FluidHeatConduction
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, return the effective heat conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     */
    template <class FluidState>
    static Scalar heatConductivity(const Params &params,
                                   const FluidState &fluidState)
    {
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);
        return FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);
    }
};
}

#endif
