// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2008-2012 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \brief Test for the non-isothermal compositional model based on flash
 *        calculations.
 */
#include "config.h"

#include <ewoms/common/quad.hh>
#include <ewoms/common/start.hh>
#include <ewoms/models/flash/flashmodel.hh>
#include "problems/co2injectionflash.hh"
#include "problems/co2injectionproblem.hh"

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(Co2InjectionFlashNIProblem,
             INHERITS_FROM(VcfvFlash, Co2InjectionBaseProblem));

SET_BOOL_PROP(Co2InjectionFlashNIProblem, EnableEnergy, true);

// for the flash model we want to use thermodynamic hints or it will
// get _very_ slow.
SET_BOOL_PROP(Co2InjectionFlashNIProblem, EnableHints, true);

// use the CO2 injection problem adapted flash solver
SET_TYPE_PROP(
    Co2InjectionFlashNIProblem, FlashSolver,
    Ewoms::Co2InjectionFlash<typename GET_PROP_TYPE(TypeTag, Scalar),
                             typename GET_PROP_TYPE(TypeTag, FluidSystem)>);

// the flash model has serious problems with the numerical
// precision. if quadruple precision math is available, we use it,
// else we increase the tolerance of the Newton solver
#if HAVE_QUAD
SET_TYPE_PROP(Co2InjectionFlashNIProblem, Scalar, quad);
#else
SET_SCALAR_PROP(Co2InjectionFlashNIProblem, NewtonRelativeTolerance, 1e-5);
#endif
}
}

int main(int argc, char **argv)
{
    typedef TTAG(Co2InjectionFlashNIProblem) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
