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
 * \brief Two-phase test for the immiscible model which uses the element-centered finite
 *        volume discretization in conjunction with automatic differentiation
 */
#ifndef EWOMS_LENS_IMMISCIBLE_ECFV_AD_HH
#define EWOMS_LENS_IMMISCIBLE_ECFV_AD_HH

#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include "problems/lensproblem.hh"

namespace Opm::Properties {

// Create new type tags
namespace TTag {
struct LensProblemEcfvAd { using InheritsFrom = std::tuple<LensBaseProblem, ImmiscibleTwoPhaseModel>; };
} // end namespace TTag

// use the element centered finite volume spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::LensProblemEcfvAd> { using type = TTag::EcfvDiscretization; };

// use automatic differentiation for this simulator
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::LensProblemEcfvAd> { using type = TTag::AutoDiffLocalLinearizer; };

// this problem works fine if the linear solver uses single precision scalars
template<class TypeTag>
struct LinearSolverScalar<TypeTag, TTag::LensProblemEcfvAd> { using type = float; };

} // namespace Opm::Properties

#endif // EWOMS_LENS_IMMISCIBLE_ECFV_AD_HH
