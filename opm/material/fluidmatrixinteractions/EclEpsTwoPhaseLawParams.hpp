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
 * \copydoc Opm::EclEpsTwoPhaseLawParams
 */
#ifndef OPM_ECL_EPS_TWO_PHASE_LAW_PARAMS_HPP
#define OPM_ECL_EPS_TWO_PHASE_LAW_PARAMS_HPP

#include "EclEpsConfig.hpp"
#include "EclEpsScalingPoints.hpp"

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

#include <string>
#include <memory>
#include <cassert>
#include <algorithm>

#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief A default implementation of the parameters for the material law adapter class
 *        which implements ECL endpoint scaleing .
 */
template <class EffLawT>
class EclEpsTwoPhaseLawParams : public EnsureFinalized
{
    typedef typename EffLawT::Params EffLawParams;
    typedef typename EffLawParams::Traits::Scalar Scalar;

public:
    typedef typename EffLawParams::Traits Traits;
    typedef Opm::EclEpsScalingPoints<Scalar> ScalingPoints;

    EclEpsTwoPhaseLawParams()
    {
    }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
#ifndef NDEBUG
        assert(config_);
        if (config_->enableSatScaling()) {
            assert(unscaledPoints_);
        }
        assert(effectiveLawParams_);
#endif
        EnsureFinalized :: finalize();
    }

    /*!
     * \brief Set the endpoint scaling configuration object.
     */
    void setConfig(std::shared_ptr<EclEpsConfig> value)
    { config_ = value; }

    /*!
     * \brief Returns the endpoint scaling configuration object.
     */
    const EclEpsConfig& config() const
    { return *config_; }

    /*!
     * \brief Set the scaling points which are seen by the nested material law
     */
    void setUnscaledPoints(std::shared_ptr<ScalingPoints> value)
    { unscaledPoints_ = value; }

    /*!
     * \brief Returns the scaling points which are seen by the nested material law
     */
    const ScalingPoints& unscaledPoints() const
    { return *unscaledPoints_; }

    /*!
     * \brief Set the scaling points which are seen by the physical model
     */
    void setScaledPoints(std::shared_ptr<ScalingPoints> value)
    { scaledPoints_ = *value; }

    /*!
     * \brief Returns the scaling points which are seen by the physical model
     */
    const ScalingPoints& scaledPoints() const
    { return scaledPoints_; }

    /*!
     * \brief Returns the scaling points which are seen by the physical model
     */
    ScalingPoints& scaledPoints()
    { return scaledPoints_; }

    /*!
     * \brief Sets the parameter object for the effective/nested material law.
     */
    void setEffectiveLawParams(std::shared_ptr<EffLawParams> value)
    { effectiveLawParams_ = value; }

    /*!
     * \brief Returns the parameter object for the effective/nested material law.
     */
    const EffLawParams& effectiveLawParams() const
    { return *effectiveLawParams_; }

private:
    std::shared_ptr<EffLawParams> effectiveLawParams_;

    std::shared_ptr<EclEpsConfig> config_;
    std::shared_ptr<ScalingPoints> unscaledPoints_;
    ScalingPoints scaledPoints_;
};

} // namespace Opm

#endif
