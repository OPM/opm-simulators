// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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
*/
/*!
 * \file
 * \copydoc Opm::EclEpsConfig
 */
#ifndef OPM_ECL_EPS_CONFIG_HPP
#define OPM_ECL_EPS_CONFIG_HPP

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

#include <string>
#include <cassert>
#include <algorithm>

namespace Opm {
/*!
 * \brief Specified which fluids are involved in a given twophase material law for
 *        endpoint scaling.
 */
enum EclTwoPhaseSystemType {
    EclGasOilSystem,
    EclOilWaterSystem,
};

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specifies the configuration used by the endpoint scaling code
 *
 * This means which quantities are scaled and how this is to be done.
 */
class EclEpsConfig
{
public:
    EclEpsConfig()
    {
        // by default, do not scale anything
        enableSatScaling_ = false;
        enablePcScaling_ = false;
        enableKrwScaling_ = false;
        enableKrnScaling_ = false;
        enableThreePointKrSatScaling_ = false;
    }

    /*!
     * \brief Specify whether saturation scaling is enabled.
     */
    void setEnableSatScaling(bool yesno)
    { enableSatScaling_ = yesno; }

    /*!
     * \brief Returns whether saturation scaling is enabled.
     */
    bool enableSatScaling() const
    { return enableSatScaling_; }

    /*!
     * \brief Specify whether three point saturation scaling is enabled for the relative
     *        permeabilities.
     */
    void setEnableThreePointKrSatScaling(bool yesno)
    { enableThreePointKrSatScaling_ = yesno; }

    /*!
     * \brief Returns whether three point saturation scaling is enabled for the relative
     *        permeabilities.
     */
    bool enableThreePointKrSatScaling() const
    { return enableThreePointKrSatScaling_; }

    /*!
     * \brief Specify whether relative permeability scaling is enabled for the wetting phase.
     */
    void setEnableKrwScaling(bool yesno)
    { enableKrwScaling_ = yesno; }

    /*!
     * \brief Returns whether relative permeability scaling is enabled for the wetting phase.
     */
    bool enableKrwScaling() const
    { return enableKrwScaling_; }

    /*!
     * \brief Specify whether relative permeability scaling is enabled for the non-wetting phase.
     */
    void setEnableKrnScaling(bool yesno)
    { enableKrnScaling_ = yesno; }

    /*!
     * \brief Returns whether relative permeability scaling is enabled for the non-wetting phase.
     */
    bool enableKrnScaling() const
    { return enableKrnScaling_; }

    /*!
     * \brief Specify whether capillary pressure scaling is enabled.
     */
    void setEnablePcScaling(bool yesno)
    { enablePcScaling_ = yesno; }

    /*!
     * \brief Returns whether capillary pressure scaling is enabled.
     */
    bool enablePcScaling() const
    { return enablePcScaling_; }

#if HAVE_OPM_PARSER
    /*!
     * \brief Reads all relevant material parameters form a cell of a parsed ECL deck.
     *
     * This requires that the opm-parser module is available.
     */
    void initFromDeck(Opm::DeckConstPtr deck,
                      Opm::EclipseStateConstPtr eclState,
                      Opm::EclTwoPhaseSystemType twoPhaseSystemType)
    {
        // find out if endpoint scaling is used in the first place
        if (!deck->hasKeyword("ENDSCALE")) {
            // it is not used, i.e., just set all enable$Foo attributes to 0 and be done
            // with it.
            enableSatScaling_ = false;
            enableThreePointKrSatScaling_ = false;
            enablePcScaling_ = false;
            enableKrwScaling_ = false;
            enableKrnScaling_ = false;
            return;
        }

        // endpoint scaling is used, i.e., at least saturation scaling needs to be enabled
        enableSatScaling_ = true;

        // check if three-point scaling is to be used for the saturations of the relative
        // permeabilities
        if (deck->hasKeyword("SCALECRS")) {
            // if the deck features the SCALECRS keyword, it must be set to 'YES'
            Opm::DeckKeywordConstPtr scalecrsKeyword = deck->getKeyword("SCALECRS");
            std::string scalecrsValue =
                scalecrsKeyword->getRecord(0)->getItem("VALUE")->getString(0);
            // convert the value of the SCALECRS keyword to upper case, just to be sure
            std::transform(scalecrsValue.begin(),
                           scalecrsValue.end(),
                           scalecrsValue.begin(),
                           ::toupper);

            if (scalecrsValue == "YES" || scalecrsValue == "Y")
                enableThreePointKrSatScaling_ = true;
            else
                enableThreePointKrSatScaling_ = false;
        }
        else
            enableThreePointKrSatScaling_ = false;

        // check if we are supposed to scale the Y axis of the capillary pressure
        if (twoPhaseSystemType == EclOilWaterSystem)
            enablePcScaling_ =
                eclState->hasDoubleGridProperty("PCW")
                || eclState->hasDoubleGridProperty("SWATINIT");

        else {
            assert(twoPhaseSystemType == EclGasOilSystem);
            enablePcScaling_ = eclState->hasDoubleGridProperty("PCG");
        }

        // check if we are supposed to scale the Y axis of the wetting phase relperm
        if (twoPhaseSystemType == EclOilWaterSystem)
            enableKrwScaling_ = eclState->hasDoubleGridProperty("KRW");
        else {
            assert(twoPhaseSystemType == EclGasOilSystem);
            enableKrwScaling_ = eclState->hasDoubleGridProperty("KRO");
        }

        // check if we are supposed to scale the Y axis of the non-wetting phase relperm
        if (twoPhaseSystemType == EclOilWaterSystem)
            enableKrnScaling_ = eclState->hasDoubleGridProperty("KRO");
        else {
            assert(twoPhaseSystemType == EclGasOilSystem);
            enableKrnScaling_ = eclState->hasDoubleGridProperty("KRG");
        }
    }
#endif

private:
    // enable scaling of the input saturations (i.e., rescale the x-Axis)
    bool enableSatScaling_;

    // use three (instead of two) points to scale the saturations for the relative
    // permeabilities.
    //
    // This means that two piecewise linear functions are used for saturation scaling
    // instead of a single linear one
    bool enableThreePointKrSatScaling_;

    // enable the scaling of the capillary pressure and relative permeability outputs
    // (i.e., rescale the y-Axis)
    bool enablePcScaling_;
    bool enableKrwScaling_;
    bool enableKrnScaling_;
};

} // namespace Opm

#endif
