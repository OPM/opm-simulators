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
 * \copydoc Opm::EclEpsTwoPhaseLawPoints
 */
#ifndef OPM_ECL_EPS_GRID_PROPERTIES_HPP
#define OPM_ECL_EPS_GRID_PROPERTIES_HPP

#include "EclEpsConfig.hpp"

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/GridProperty.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgfnTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SlgofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Sof3Table.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SwfnTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SwofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

#include <opm/material/common/Means.hpp>

#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>

namespace Opm {
/*!
 * \brief Collects all grid properties which are relevant for end point scaling.
 *
 * This class is used for both, the drainage and the imbibition variants of the ECL
 * keywords.
 */
class EclEpsGridProperties
{
    typedef std::vector<int> IntData;
    typedef std::vector<double> DoubleData;

public:
#if HAVE_ECL_INPUT
    void initFromDeck(const Opm::Deck& /* deck */,
                      const Opm::EclipseState& eclState,
                      bool useImbibition)
    {
        std::string kwPrefix = useImbibition?"I":"";

        if (useImbibition)
            satnum = &eclState.get3DProperties().getIntGridProperty("IMBNUM").getData();
        else
            satnum = &eclState.get3DProperties().getIntGridProperty("SATNUM").getData();

        retrieveGridPropertyData_(&swl, eclState, kwPrefix+"SWL");
        retrieveGridPropertyData_(&sgl, eclState, kwPrefix+"SGL");
        retrieveGridPropertyData_(&swcr, eclState, kwPrefix+"SWCR");
        retrieveGridPropertyData_(&sgcr, eclState, kwPrefix+"SGCR");
        retrieveGridPropertyData_(&sowcr, eclState, kwPrefix+"SOWCR");
        retrieveGridPropertyData_(&sogcr, eclState, kwPrefix+"SOGCR");
        retrieveGridPropertyData_(&swu, eclState, kwPrefix+"SWU");
        retrieveGridPropertyData_(&sgu, eclState, kwPrefix+"SGU");
        retrieveGridPropertyData_(&pcw, eclState, kwPrefix+"PCW");
        retrieveGridPropertyData_(&pcg, eclState, kwPrefix+"PCG");
        retrieveGridPropertyData_(&krw, eclState, kwPrefix+"KRW");
        retrieveGridPropertyData_(&kro, eclState, kwPrefix+"KRO");
        retrieveGridPropertyData_(&krg, eclState, kwPrefix+"KRG");

        // _may_ be needed to calculate the Leverett capillary pressure scaling factor
        const auto& ecl3dProps = eclState.get3DProperties();
        poro = &ecl3dProps.getDoubleGridProperty("PORO").getData();

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMX")) {
            permx = &ecl3dProps.getDoubleGridProperty("PERMX").getData();
            permy = permx;
            permz = permx;
        }

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMY"))
            permy = &ecl3dProps.getDoubleGridProperty("PERMY").getData();

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMZ"))
            permz = &ecl3dProps.getDoubleGridProperty("PERMZ").getData();
    }
#endif

    const IntData* satnum;

    const DoubleData* swl;
    const DoubleData* sgl;
    const DoubleData* swcr;
    const DoubleData* sgcr;
    const DoubleData* sowcr;
    const DoubleData* sogcr;
    const DoubleData* swu;
    const DoubleData* sgu;
    const DoubleData* pcw;
    const DoubleData* pcg;
    const DoubleData* krw;
    const DoubleData* kro;
    const DoubleData* krg;
    const DoubleData* poro;
    const DoubleData* permx;
    const DoubleData* permy;
    const DoubleData* permz;

private:
#if HAVE_ECL_INPUT
    // this method makes sure that a grid property is not created if it is not explicitly
    // mentioned in the deck. (saves memory.)
    void retrieveGridPropertyData_(const DoubleData **data,
                                   const Opm::EclipseState& eclState,
                                   const std::string& properyName)
    {
        (*data) = 0;
        if (eclState.get3DProperties().hasDeckDoubleGridProperty(properyName))
            (*data) = &eclState.get3DProperties().getDoubleGridProperty(properyName).getData();
    }
#endif
};
#endif
}
