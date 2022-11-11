/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/simulators/wells/MultisegmentWellEquations.hpp>

namespace Opm {

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::clear()
{
    duneB_ = 0.0;
    duneC_ = 0.0;
    duneD_ = 0.0;
    resWell_ = 0.0;
    duneDSolver_.reset();
}

#define INSTANCE(numWellEq, numEq) \
template class MultisegmentWellEquations<double,numWellEq,numEq>;

INSTANCE(2,1)
INSTANCE(2,2)
INSTANCE(2,6)
INSTANCE(3,2)
INSTANCE(3,3)
INSTANCE(3,4)
INSTANCE(4,3)
INSTANCE(4,4)
INSTANCE(4,5)

}
