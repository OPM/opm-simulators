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

public:
#if HAVE_ECL_INPUT

    EclEpsGridProperties(const Opm::EclipseState& eclState,
                         bool useImbibition)
    {
        std::string kwPrefix = useImbibition?"I":"";
        if (useImbibition)
            global_satnum = &eclState.get3DProperties().getIntGridProperty("IMBNUM").getData();
        else
            global_satnum = &eclState.get3DProperties().getIntGridProperty("SATNUM").getData();

        retrieveGridPropertyData_(&this->global_swl, eclState, kwPrefix+"SWL");
        retrieveGridPropertyData_(&this->global_sgl, eclState, kwPrefix+"SGL");
        retrieveGridPropertyData_(&this->global_swcr, eclState, kwPrefix+"SWCR");
        retrieveGridPropertyData_(&this->global_sgcr, eclState, kwPrefix+"SGCR");
        retrieveGridPropertyData_(&this->global_sowcr, eclState, kwPrefix+"SOWCR");
        retrieveGridPropertyData_(&this->global_sogcr, eclState, kwPrefix+"SOGCR");
        retrieveGridPropertyData_(&this->global_swu, eclState, kwPrefix+"SWU");
        retrieveGridPropertyData_(&this->global_sgu, eclState, kwPrefix+"SGU");
        retrieveGridPropertyData_(&this->global_pcw, eclState, kwPrefix+"PCW");
        retrieveGridPropertyData_(&this->global_pcg, eclState, kwPrefix+"PCG");
        retrieveGridPropertyData_(&this->global_krw, eclState, kwPrefix+"KRW");
        retrieveGridPropertyData_(&this->global_kro, eclState, kwPrefix+"KRO");
        retrieveGridPropertyData_(&this->global_krg, eclState, kwPrefix+"KRG");

        // _may_ be needed to calculate the Leverett capillary pressure scaling factor
        const auto& ecl3dProps = eclState.get3DProperties();
        this->global_poro = &ecl3dProps.getDoubleGridProperty("PORO").getData();

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMX")) {
            this->global_permx = &ecl3dProps.getDoubleGridProperty("PERMX").getData();
            this->global_permy = this->global_permx;
            this->global_permz = this->global_permx;
        }

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMY"))
            this->global_permy = &ecl3dProps.getDoubleGridProperty("PERMY").getData();

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMZ"))
            this->global_permz = &ecl3dProps.getDoubleGridProperty("PERMZ").getData();
    }
#endif



    unsigned satRegion(std::size_t global_index) const {
        return this->global_satnum->operator[](global_index) - 1;
    }

    double permx(std::size_t global_index) const {
        return this->global_permx->operator[](global_index);
    }

    double permy(std::size_t global_index) const {
        return this->global_permy->operator[](global_index);
    }

    double permz(std::size_t global_index) const {
        return this->global_permz->operator[](global_index);
    }

    double poro(std::size_t global_index) const {
        return this->global_poro->operator[](global_index);
    }

    const double * swl(std::size_t global_index) const {
        return this->satfunc(this->global_swl, global_index);
    }

    const double * sgl(std::size_t global_index) const {
        return this->satfunc(this->global_sgl, global_index);
    }

    const double * swcr(std::size_t global_index) const {
        return this->satfunc(this->global_swcr, global_index);
    }

    const double * sgcr(std::size_t global_index) const {
        return this->satfunc(this->global_sgcr, global_index);
    }

    const double * sowcr(std::size_t global_index) const {
        return this->satfunc(this->global_sowcr, global_index);
    }

    const double * sogcr(std::size_t global_index) const {
        return this->satfunc(this->global_sogcr, global_index);
    }

    const double * swu(std::size_t global_index) const {
        return this->satfunc(this->global_swu, global_index);
    }

    const double * sgu(std::size_t global_index) const {
        return this->satfunc(this->global_sgu, global_index);
    }

    const double * pcw(std::size_t global_index) const {
        return this->satfunc(this->global_pcw, global_index);
    }

    const double * pcg(std::size_t global_index) const {
        return this->satfunc(this->global_pcg, global_index);
    }

    const double * krw(std::size_t global_index) const {
        return this->satfunc(this->global_krw, global_index);
    }

    const double * krg(std::size_t global_index) const {
        return this->satfunc(this->global_krg, global_index);
    }

    const double * kro(std::size_t global_index) const {
        return this->satfunc(this->global_kro, global_index);
    }


private:
#if HAVE_ECL_INPUT
    // this method makes sure that a grid property is not created if it is not explicitly
    // mentioned in the deck. (saves memory.)
    void retrieveGridPropertyData_(const std::vector<double> **data,
                                   const Opm::EclipseState& eclState,
                                   const std::string& properyName)
    {
        (*data) = 0;
        if (eclState.get3DProperties().hasDeckDoubleGridProperty(properyName))
            (*data) = &eclState.get3DProperties().getDoubleGridProperty(properyName).getData();
    }
#endif

    const double * satfunc(const std::vector<double>* data, std::size_t global_index) const {
        if (data)
            return &(data->operator[](global_index));
        return nullptr;
    }


    const std::vector<int>*    global_satnum;
    const std::vector<double>* global_swl;
    const std::vector<double>* global_sgl;
    const std::vector<double>* global_swcr;
    const std::vector<double>* global_sgcr;
    const std::vector<double>* global_sowcr;
    const std::vector<double>* global_sogcr;
    const std::vector<double>* global_swu;
    const std::vector<double>* global_sgu;
    const std::vector<double>* global_pcw;
    const std::vector<double>* global_pcg;
    const std::vector<double>* global_krw;
    const std::vector<double>* global_kro;
    const std::vector<double>* global_krg;

    const std::vector<double>* global_permx;
    const std::vector<double>* global_permy;
    const std::vector<double>* global_permz;
    const std::vector<double>* global_poro;
};
}
#endif

