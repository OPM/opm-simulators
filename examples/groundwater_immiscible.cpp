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
 * \brief Test for the immisicible VCVF discretization with only a single phase
 */
#include "config.h"

#include <opm/models/common/darcyfluxmodule.hh>
#include <opm/models/discretization/common/fvbasefdlocallinearizer.hh>
#include <opm/models/discretization/vcfv/vcfvdiscretization.hh>
#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/models/io/dgfvanguard.hh>
#include <opm/models/utils/start.hh>

#include "problems/groundwaterproblem.hh"

namespace Opm::Properties {

// Create new type tags
namespace TTag {

struct GroundWaterProblem
{ using InheritsFrom = std::tuple<GroundWaterBaseProblem, ImmiscibleSinglePhaseModel>; };

} // end namespace TTag

//! We use a vertex centered finite volume method
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::GroundWaterProblem>
{ using type = TTag::VcfvDiscretization; };

//! Use the Darcy relation to determine the phase velocity
template<class TypeTag>
struct FluxModule<TypeTag, TTag::GroundWaterProblem>
{ using type = DarcyFluxModule<TypeTag>; };

// //! Use finite differences to linearize the system of PDEs
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::GroundWaterProblem>
{ using type = TTag::FiniteDifferenceLocalLinearizer; };

} // namespace Opm::Properties

int main(int argc, char **argv)
{
    using ProblemTypeTag = Opm::Properties::TTag::GroundWaterProblem;
    return Opm::start<ProblemTypeTag>(argc, argv, true);
}
