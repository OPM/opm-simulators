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

/*
 * \file
 *
 * \brief Test for the isothermal primary variable switching model
 *        using "obstacle" problem and the VCVF discretization.
 */
#include "config.h"

#include <ewoms/common/start.hh>
#include <ewoms/models/pvs/pvsmodel.hh>
#include "problems/obstacleproblem.hh"

BEGIN_PROPERTIES

NEW_TYPE_TAG(ObstacleProblem, INHERITS_FROM(PvsModel, ObstacleBaseProblem));

// Verbosity of the PVS model (0=silent, 1=medium, 2=chatty)
SET_INT_PROP(ObstacleProblem, PvsVerbosity, 1);

END_PROPERTIES

int main(int argc, char **argv)
{
    typedef TTAG(ObstacleProblem) ProblemTypeTag;
    return Opm::start<ProblemTypeTag>(argc, argv);
}
