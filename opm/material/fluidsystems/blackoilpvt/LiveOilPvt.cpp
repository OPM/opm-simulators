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

#include "config.h"

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt_impl.hpp>

namespace Opm {

template class LiveOilPvt<double>;

#define INSTANCE_INT_ENERGY(...) \
    template __VA_ARGS__ \
    LiveOilPvt<double>::internalEnergy<__VA_ARGS__>(unsigned, \
                                                    const __VA_ARGS__&, \
                                                    const __VA_ARGS__&, \
                                                    const __VA_ARGS__&) const;

INSTANCE_INT_ENERGY(DenseAd::Evaluation<double,1,0u>)
INSTANCE_INT_ENERGY(DenseAd::Evaluation<double,2,0u>)
INSTANCE_INT_ENERGY(DenseAd::Evaluation<double,3,0u>)
INSTANCE_INT_ENERGY(DenseAd::Evaluation<double,4,0u>)
INSTANCE_INT_ENERGY(DenseAd::Evaluation<double,5,0u>)
INSTANCE_INT_ENERGY(DenseAd::Evaluation<double,6,0u>)

#define INSTANCE_SAT_GAS(...) \
    template __VA_ARGS__ \
    LiveOilPvt<double>::saturatedGasDissolutionFactor<__VA_ARGS__>(unsigned, \
                                                                   const __VA_ARGS__&, \
                                                                   const __VA_ARGS__&, \
                                                                   const __VA_ARGS__&, \
                                                                   __VA_ARGS__) const; \
    template __VA_ARGS__ \
    LiveOilPvt<double>::saturationPressure<__VA_ARGS__>(unsigned, \
                                                        const __VA_ARGS__&, \
                                                        const __VA_ARGS__&) const; \
    template __VA_ARGS__ \
    LiveOilPvt<double>::diffusionCoefficient<__VA_ARGS__>(const __VA_ARGS__&, \
                                                          const __VA_ARGS__&, \
                                                          unsigned) const;

INSTANCE_SAT_GAS(double)
INSTANCE_SAT_GAS(DenseAd::Evaluation<double,1,0u>)
INSTANCE_SAT_GAS(DenseAd::Evaluation<double,2,0u>)
INSTANCE_SAT_GAS(DenseAd::Evaluation<double,3,0u>)
INSTANCE_SAT_GAS(DenseAd::Evaluation<double,4,0u>)
INSTANCE_SAT_GAS(DenseAd::Evaluation<double,5,0u>)
INSTANCE_SAT_GAS(DenseAd::Evaluation<double,6,0u>)

#define INSTANCE_FULL(...) \
    template __VA_ARGS__ \
    LiveOilPvt<double>::viscosity<__VA_ARGS__>(unsigned, \
                                               const __VA_ARGS__&, \
                                               const __VA_ARGS__&, \
                                               const __VA_ARGS__&) const; \
    template __VA_ARGS__ \
    LiveOilPvt<double>::saturatedViscosity<__VA_ARGS__>(unsigned, \
                                                        const __VA_ARGS__&, \
                                                        const __VA_ARGS__&) const; \
    template __VA_ARGS__ \
    LiveOilPvt<double>::inverseFormationVolumeFactor<__VA_ARGS__>(unsigned, \
                                                                  const __VA_ARGS__&, \
                                                                  const __VA_ARGS__&, \
                                                                  const __VA_ARGS__&) const; \
    template __VA_ARGS__ \
    LiveOilPvt<double>::saturatedInverseFormationVolumeFactor<__VA_ARGS__>(unsigned, \
                                                                           const __VA_ARGS__&, \
                                                                           const __VA_ARGS__&) const; \
    template __VA_ARGS__ \
    LiveOilPvt<double>::saturatedGasDissolutionFactor<__VA_ARGS__>(unsigned, \
                                                                   const __VA_ARGS__&, \
                                                                   const __VA_ARGS__&) const;

INSTANCE_FULL(double)
INSTANCE_FULL(DenseAd::Evaluation<double,1,0u>)
INSTANCE_FULL(DenseAd::Evaluation<double,2,0u>)
INSTANCE_FULL(DenseAd::Evaluation<double,3,0u>)
INSTANCE_FULL(DenseAd::Evaluation<double,4,0u>)
INSTANCE_FULL(DenseAd::Evaluation<double,5,0u>)
INSTANCE_FULL(DenseAd::Evaluation<double,6,0u>)
INSTANCE_FULL(DenseAd::Evaluation<double,7,0u>)
INSTANCE_FULL(DenseAd::Evaluation<double,8,0u>)
INSTANCE_FULL(DenseAd::Evaluation<double,9,0u>)

} // namespace Opm
