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
 * \brief Test for the non-isothermal compositional model based on flash
 *        calculations.
 */
#include "config.h"

#include <opm/material/common/quad.hpp>
#include <opm/models/utils/start.hh>
#include <opm/models/flash/flashmodel.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include "problems/co2injectionflash.hh"
#include "problems/co2injectionproblem.hh"

namespace Opm::Properties {

// Create new type tags
namespace TTag {
struct Co2InjectionFlashNiEcfvProblem { using InheritsFrom = std::tuple<Co2InjectionBaseProblem, FlashModel>; };
} // end namespace TTag
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::Co2InjectionFlashNiEcfvProblem> { using type = TTag::EcfvDiscretization; };

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::Co2InjectionFlashNiEcfvProblem> { static constexpr bool value = true; };

//! Use automatic differentiation to linearize the system of PDEs
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::Co2InjectionFlashNiEcfvProblem> { using type = TTag::AutoDiffLocalLinearizer; };

// use the CO2 injection problem adapted flash solver
template<class TypeTag>
struct FlashSolver<TypeTag, TTag::Co2InjectionFlashNiEcfvProblem>
{ using type = Opm::Co2InjectionFlash<GetPropType<TypeTag, Properties::Scalar>,
                                      GetPropType<TypeTag, Properties::FluidSystem>>; };

// the flash model has serious problems with the numerical
// precision. if quadruple precision math is available, we use it,
// else we increase the tolerance of the Newton solver
#if HAVE_QUAD
template<class TypeTag>
struct Scalar<TypeTag, TTag::Co2InjectionFlashNiEcfvProblem> { using type = quad; };

// the default linear solver used for this problem (-> AMG) cannot be used with quadruple
// precision scalars... (this seems to only apply to Dune >= 2.4)
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::Co2InjectionFlashNiEcfvProblem> { using type = TTag::ParallelBiCGStabLinearSolver; };
#else
template<class TypeTag>
struct NewtonTolerance<TypeTag, TTag::Co2InjectionFlashNiEcfvProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-5;
};
#endif

} // namespace Opm::Properties

int main(int argc, char **argv)
{
    using EcfvProblemTypeTag = Opm::Properties::TTag::Co2InjectionFlashNiEcfvProblem;
    return Opm::start<EcfvProblemTypeTag>(argc, argv);
}
