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
 * \brief Simulation of the injection problem using the VCVF discretization
 *        assuming immisicibility and with energy enabled.
 */
#include "config.h"

#include <opm/models/utils/start.hh>
#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/models/discretization/vcfv/vcfvdiscretization.hh>
#include "problems/co2injectionproblem.hh"

namespace Opm::Properties {

namespace TTag {

struct Co2InjectionImmiscibleNiVcfvProblem
{ using InheritsFrom = std::tuple<Co2InjectionBaseProblem, ImmiscibleModel>; };

} // namespace TTag

template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::Co2InjectionImmiscibleNiVcfvProblem> { using type = TTag::VcfvDiscretization; };

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::Co2InjectionImmiscibleNiVcfvProblem> { static constexpr bool value = true; };

} // namespace Opm::Properties

////////////////////////
// the main function
////////////////////////
int main(int argc, char **argv)
{
    using VcfvProblemTypeTag = Opm::Properties::TTag::Co2InjectionImmiscibleNiVcfvProblem;
    return Opm::start<VcfvProblemTypeTag>(argc, argv);
}
