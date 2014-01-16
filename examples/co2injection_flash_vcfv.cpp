/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
 * \brief Test for the isothermal compositional model based on flash
 *        calculations.
 */
#include "config.h"

#if HAVE_QUAD
#include <ewoms/common/quad.hh>
#endif

#include <ewoms/common/start.hh>
#include <ewoms/models/flash/flashmodel.hh>
#include <ewoms/disc/vcfv/vcfvdiscretization.hh>
#include "problems/co2injectionflash.hh"
#include "problems/co2injectionproblem.hh"

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(Co2InjectionFlashVcfvProblem, INHERITS_FROM(FlashModel, Co2InjectionBaseProblem));
SET_TAG_PROP(Co2InjectionFlashVcfvProblem, SpatialDiscretizationSplice, VcfvDiscretization);

// for the flash model we want to use thermodynamic hints or it will
// get _very_ slow.
SET_BOOL_PROP(Co2InjectionFlashVcfvProblem, EnableHints, true);

// use the flash solver adapted to the CO2 injection problem
SET_TYPE_PROP(
    Co2InjectionFlashVcfvProblem, FlashSolver,
    Ewoms::Co2InjectionFlash<typename GET_PROP_TYPE(TypeTag, Scalar),
                             typename GET_PROP_TYPE(TypeTag, FluidSystem)>);

// the flash model has serious problems with the numerical
// precision. if quadruple precision math is available, we use it,
// else we increase the tolerance of the Newton solver
#if HAVE_QUAD
SET_TYPE_PROP(Co2InjectionFlashVcfvProblem, Scalar, quad);
#else
SET_SCALAR_PROP(Co2InjectionFlashVcfvProblem, NewtonRelativeTolerance, 1e-5);
#endif
} // namespace Properties
} // namespace Opm

int main(int argc, char **argv)
{
    typedef TTAG(Co2InjectionFlashVcfvProblem) VcfvProblemTypeTag;
    return Ewoms::start<VcfvProblemTypeTag>(argc, argv);
}
