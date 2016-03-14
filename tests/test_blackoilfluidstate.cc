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
 * \brief This test makes sure that the programming interface is
 *        observed by all fluid systems
 */
#include "config.h"

#include <opm/material/localad/Evaluation.hpp>
#include <opm/material/localad/Math.hpp>
#include <opm/material/checkFluidSystem.hpp>

#include <ewoms/models/blackoil/blackoilmodel.hh>
#include <ewoms/models/blackoil/blackoilfluidstate.hh>

#include <opm/common/utility/platform_dependent/disable_warnings.h>

// include dune's MPI helper header
#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

// check the API of the black-oil fluid states
template <class Scalar>
void testBlackOilFluidState()
{
    typedef TTAG(BlackOilModel) TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef Ewoms::BlackOilFluidState<TypeTag> FluidState;

    FluidState fs;
    checkFluidState<Evaluation>(fs);
}

class TestAdTag;

template <class Scalar>
inline void testAll()
{
    typedef Opm::LocalAd::Evaluation<Scalar, TestAdTag, 3> Evaluation;

    // ensure that all fluid states are API-compliant
    testBlackOilFluidState<Scalar>();
    testBlackOilFluidState<Evaluation>();
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);
    testAll<double>();
    testAll<float>()
;
    return 0;
}
