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

#include <config.h>
#include <opm/models/blackoil/blackoilmicpparams.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/MICPpara.hpp>
#endif

#include <algorithm>
#include <stdexcept>
#include <type_traits>

namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
template<bool enableMICP>
void BlackOilMICPParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // some sanity checks: if MICP is enabled, the MICP keyword must be
    // present, if MICP is disabled the keyword must not be present.
    if constexpr (enableMICP) {
        if (!eclState.runspec().micp()) {
            throw std::runtime_error("Non-trivial MICP treatment requested at compile time, but "
                                     "the deck does not contain the MICP keyword");
        }
    }
    else {
        if (eclState.runspec().micp()) {
            throw std::runtime_error("MICP treatment disabled at compile time, but the deck "
                                     "contains the MICP keyword");
        }
    }

    if (!eclState.runspec().micp())
        return; // MICP treatment is supposed to be disabled*/

    // initialize the objects which deal with the MICPpara keyword
    const auto& MICPpara = eclState.getMICPpara();
    densityBiofilm_ = MICPpara.getDensityBiofilm();
    densityCalcite_ = MICPpara.getDensityCalcite();
    detachmentRate_ = MICPpara.getDetachmentRate();
    criticalPorosity_ = MICPpara.getCriticalPorosity();
    fittingFactor_ = MICPpara.getFittingFactor();
    halfVelocityOxygen_ = MICPpara.getHalfVelocityOxygen();
    halfVelocityUrea_ = MICPpara.getHalfVelocityUrea();
    maximumGrowthRate_ = MICPpara.getMaximumGrowthRate();
    maximumUreaUtilization_ = MICPpara.getMaximumUreaUtilization();
    microbialAttachmentRate_ = MICPpara.getMicrobialAttachmentRate();
    microbialDeathRate_ = MICPpara.getMicrobialDeathRate();
    minimumPermeability_ = MICPpara.getMinimumPermeability();
    oxygenConsumptionFactor_ = MICPpara.getOxygenConsumptionFactor();
    yieldGrowthCoefficient_ = MICPpara.getYieldGrowthCoefficient();
    maximumOxygenConcentration_ = MICPpara.getMaximumOxygenConcentration();
    maximumUreaConcentration_ = MICPpara.getMaximumUreaConcentration();
    toleranceBeforeClogging_ = MICPpara.getToleranceBeforeClogging();

    // obtain the porosity for the clamp in the blackoilnewtonmethod
    if constexpr (std::is_same_v<Scalar, float>) {
        const auto phi = eclState.fieldProps().get_double("PORO");
        phi_.resize(phi.size());
        std::copy(phi.begin(), phi.end(), phi_.begin());
    } else {
        phi_ = eclState.fieldProps().get_double("PORO");
    }
}
#endif

#define INSTANTIATE_TYPE(T)                                                         \
    template struct BlackOilMICPParams<T>;                                          \
    template void BlackOilMICPParams<T>::initFromState<false>(const EclipseState&); \
    template void BlackOilMICPParams<T>::initFromState<true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
