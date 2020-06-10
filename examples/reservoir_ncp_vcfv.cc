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
 * \brief Test for the reservoir problem using the NCP model, the VCFV discretization and
 *        finite differences.
 */
#include "config.h"

#include <opm/models/utils/start.hh>
#include <opm/models/ncp/ncpmodel.hh>
#include <opm/models/discretization/vcfv/vcfvdiscretization.hh>
#include "problems/reservoirproblem.hh"

namespace Opm::Properties {

// Create new type tags
namespace TTag {
struct ReservoirNcpVcfvProblem { using InheritsFrom = std::tuple<ReservoirBaseProblem, NcpModel>; };
} // end namespace TTag

// Select the vertex centered finite volume method as spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::ReservoirNcpVcfvProblem> { using type = TTag::VcfvDiscretization; };

// enable the storage cache for this problem so that the storage cache receives wider
// testing
template<class TypeTag>
struct EnableStorageCache<TypeTag, TTag::ReservoirNcpVcfvProblem> { static constexpr bool value = true; };

// reduce the base epsilon for the finite difference method to 10^-11. for some reason
// the simulator converges better with this. (TODO: use automatic differentiation?)
template<class TypeTag>
struct BaseEpsilon<TypeTag, TTag::ReservoirNcpVcfvProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-11;
};


} // namespace Opm::Properties

int main(int argc, char **argv)
{
    using ProblemTypeTag = Opm::Properties::TTag::ReservoirNcpVcfvProblem;
    return Opm::start<ProblemTypeTag>(argc, argv);
}
