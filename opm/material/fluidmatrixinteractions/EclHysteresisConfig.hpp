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
 * \copydoc Opm::EclHysteresisConfig
 */
#ifndef OPM_ECL_HYSTERESIS_CONFIG_HPP
#define OPM_ECL_HYSTERESIS_CONFIG_HPP

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/Deck/DeckItem.hpp>
#endif

#include <opm/material/common/Exceptions.hpp>

#include <string>
#include <cassert>
#include <algorithm>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specifies the configuration used by the ECL kr/pC hysteresis code
 */
class EclHysteresisConfig
{
public:
    EclHysteresisConfig()
    {
        enableHysteresis_ = false;
        pcHysteresisModel_ = 0;
        krHysteresisModel_ = 0;
    }

    /*!
     * \brief Specify whether hysteresis is enabled or not.
     */
    void setEnableHysteresis(bool yesno)
    { enableHysteresis_ = yesno; }

    /*!
     * \brief Returns whether hysteresis is enabled.
     */
    bool enableHysteresis() const
    { return enableHysteresis_; }

    /*!
     * \brief Set the type of the hysteresis model which is used for capillary pressure.
     *
     * -1: capillary pressure hysteresis is disabled
     * 0: use the Killough model for capillary pressure hysteresis
     */
    void setPcHysteresisModel(int value)
    { pcHysteresisModel_ = value; }

    /*!
     * \brief Return the type of the hysteresis model which is used for capillary pressure.
     *
     * -1: capillary pressure hysteresis is disabled
     * 0: use the Killough model for capillary pressure hysteresis
     */
    int pcHysteresisModel() const
    { return pcHysteresisModel_; }

    /*!
     * \brief Set the type of the hysteresis model which is used for relative permeability.
     *
     * -1: relperm hysteresis is disabled
     * 0: use the Carlson model for relative permeability hysteresis
     */
    void setKrHysteresisModel(int value)
    { krHysteresisModel_ = value; }

    /*!
     * \brief Return the type of the hysteresis model which is used for relative permeability.
     *
     * -1: relperm hysteresis is disabled
     * 0: use the Carlson model for relative permeability hysteresis
     */
    int krHysteresisModel() const
    { return krHysteresisModel_; }

#if HAVE_ECL_INPUT
    /*!
     * \brief Reads all relevant material parameters form a cell of a parsed ECL deck.
     *
     * This requires that the opm-parser module is available.
     */
    void initFromDeck(const Opm::Deck& deck)
    {
        enableHysteresis_ = false;

        if (!deck.hasKeyword("SATOPTS"))
            return;

        const auto& satoptsItem = deck.getKeyword("SATOPTS").getRecord(0).getItem(0);
        for (unsigned i = 0; i < satoptsItem.size(); ++i) {
            std::string satoptsValue = satoptsItem.get< std::string >(0);
            std::transform(satoptsValue.begin(),
                           satoptsValue.end(),
                           satoptsValue.begin(),
                           ::toupper);

            if (satoptsValue == "HYSTER")
                enableHysteresis_ = true;
        }

        // check for the (deprecated) HYST keyword
        if (deck.hasKeyword("HYST"))
            enableHysteresis_ = true;

        if (!enableHysteresis_)
            return;

        if (!deck.hasKeyword("EHYSTR"))
            throw std::runtime_error("Enabling hysteresis via the HYST parameter for SATOPTS requires the "
                                     "presence of the EHYSTR keyword");

        const auto& ehystrKeyword = deck.getKeyword("EHYSTR");
        if (deck.hasKeyword("NOHYKR"))
            krHysteresisModel_ = -1;
        else {
            krHysteresisModel_ = ehystrKeyword.getRecord(0).getItem("relative_perm_hyst").get< int >(0);
            if (krHysteresisModel_ != 0)
                throw std::runtime_error("Only the Carlson kr hystersis model (indicated by a 0 on the second item"
                                         " of the 'EHYSTR' keyword) is supported");
        }

        if (deck.hasKeyword("NOHYPC"))
            pcHysteresisModel_ = -1;
        else {
            // if capillary pressure hysteresis is enabled, Eclipse always uses the
            // Killough model
            pcHysteresisModel_ = 0;
        }
    }
#endif

private:
    // enable hysteresis at all
    bool enableHysteresis_;

    // the capillary pressure and the relperm hysteresis models to be used
    int pcHysteresisModel_;
    int krHysteresisModel_;
};

} // namespace Opm

#endif
