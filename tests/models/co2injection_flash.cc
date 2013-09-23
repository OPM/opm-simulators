// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 *
 * \brief Test for the isothermal compositional model based on flash
 *        calculations.
 */
#include "config.h"

#if HAVE_QUAD
#include <ewoms/common/quad.hh>
#endif

#include <ewoms/common/start.hh>
#include <ewoms/models/flash/flashmodel.hh>
#include "problems/co2injectionflash.hh"
#include "problems/co2injectionproblem.hh"

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(Co2InjectionFlashProblem, INHERITS_FROM(VcfvFlash, Co2InjectionBaseProblem));

// for the flash model we want to use thermodynamic hints or it will
// get _very_ slow.
SET_BOOL_PROP(Co2InjectionFlashProblem, EnableHints, true);

// use the flash solver adapted to the CO2 injection problem
SET_TYPE_PROP(Co2InjectionFlashProblem,
              FlashSolver,
              Ewoms::Co2InjectionFlash<typename GET_PROP_TYPE(TypeTag, Scalar),
              typename GET_PROP_TYPE(TypeTag, FluidSystem)> );

// the flash model has serious problems with the numerical
// precision. if quadruple precision math is available, we use it,
// else we increase the tolerance of the Newton solver
#if HAVE_QUAD
SET_TYPE_PROP(Co2InjectionFlashProblem, Scalar, quad);
#else
SET_SCALAR_PROP(Co2InjectionFlashProblem, NewtonRelativeTolerance, 1e-5);
#endif

}
}

int main(int argc, char** argv)
{
    typedef TTAG(Co2InjectionFlashProblem) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
