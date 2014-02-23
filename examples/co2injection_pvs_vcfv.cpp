/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 * \brief Test for the isothermal primary variable switching model using the VCVF
 *        discretization.
 */
#include "config.h"

#include <ewoms/common/start.hh>
#include <ewoms/models/pvs/pvsmodel.hh>
#include <ewoms/disc/vcfv/vcfvdiscretization.hh>

#include "problems/co2injectionproblem.hh"

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(Co2InjectionPvsVcfvProblem, INHERITS_FROM(PvsModel, Co2InjectionBaseProblem));
SET_TAG_PROP(Co2InjectionPvsVcfvProblem, SpatialDiscretizationSplice, VcfvDiscretization);
}}

int main(int argc, char **argv)
{
    typedef TTAG(Co2InjectionPvsVcfvProblem) VcfvProblemTypeTag;
    return Ewoms::start<VcfvProblemTypeTag>(argc, argv);
}
