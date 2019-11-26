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
#include <memory>
namespace Opm {
/*!
 * \brief Collects all grid properties which are relevant for end point scaling.
 *
 * This class is used for both, the drainage and the imbibition variants of the ECL
 * keywords.
 */

namespace {

template <typename T>
std::unique_ptr<std::vector<T>> compressed_copy(const std::vector<T>& global_vector, const std::vector<int>& compressedToCartesianElemIdx) {
    std::unique_ptr<std::vector<T>> compressed = std::make_unique<std::vector<T>>(compressedToCartesianElemIdx.size());

    for (std::size_t active_index = 0; active_index < compressedToCartesianElemIdx.size(); active_index++) {
        auto global_index = compressedToCartesianElemIdx[active_index];
        compressed->operator[](active_index) = global_vector[global_index];
    }

    return compressed;
}


std::unique_ptr<std::vector<double>> try_get(const Eclipse3DProperties& props, const std::string& keyword, const std::vector<int>& compressedToCartesianElemIdx) {
    if (props.hasDeckDoubleGridProperty(keyword))
        return compressed_copy(props.getDoubleGridProperty(keyword).getData(), compressedToCartesianElemIdx);

    return {};
}

}



class EclEpsGridProperties
{

public:
#if HAVE_ECL_INPUT

    EclEpsGridProperties(const Opm::EclipseState& eclState,
                         bool useImbibition,
                         const std::vector<int>& compressedToCartesianElemIdx)
    {
        std::string kwPrefix = useImbibition?"I":"";
        const auto& ecl3dProps = eclState.get3DProperties();

        if (useImbibition)
            compressed_satnum = compressed_copy(ecl3dProps.getIntGridProperty("IMBNUM").getData(), compressedToCartesianElemIdx);
        else
            compressed_satnum = compressed_copy(ecl3dProps.getIntGridProperty("SATNUM").getData(), compressedToCartesianElemIdx);

        this->compressed_swl = try_get( ecl3dProps, kwPrefix+"SWL", compressedToCartesianElemIdx);
        this->compressed_sgl = try_get( ecl3dProps, kwPrefix+"SGL", compressedToCartesianElemIdx);
        this->compressed_swcr = try_get( ecl3dProps, kwPrefix+"SWCR", compressedToCartesianElemIdx);
        this->compressed_sgcr = try_get( ecl3dProps, kwPrefix+"SGCR", compressedToCartesianElemIdx);
        this->compressed_sowcr = try_get( ecl3dProps, kwPrefix+"SOWCR", compressedToCartesianElemIdx);
        this->compressed_sogcr = try_get( ecl3dProps, kwPrefix+"SOGCR", compressedToCartesianElemIdx);
        this->compressed_swu = try_get( ecl3dProps, kwPrefix+"SWU", compressedToCartesianElemIdx);
        this->compressed_sgu = try_get( ecl3dProps, kwPrefix+"SGU", compressedToCartesianElemIdx);
        this->compressed_pcw = try_get( ecl3dProps, kwPrefix+"PCW", compressedToCartesianElemIdx);
        this->compressed_pcg = try_get( ecl3dProps, kwPrefix+"PCG", compressedToCartesianElemIdx);
        this->compressed_krw = try_get( ecl3dProps, kwPrefix+"KRW", compressedToCartesianElemIdx);
        this->compressed_kro = try_get( ecl3dProps, kwPrefix+"KRO", compressedToCartesianElemIdx);
        this->compressed_krg = try_get( ecl3dProps, kwPrefix+"KRG", compressedToCartesianElemIdx);

        // _may_ be needed to calculate the Leverett capillary pressure scaling factor
        if (ecl3dProps.hasDeckDoubleGridProperty("PORO"))
            this->compressed_poro = compressed_copy(ecl3dProps.getDoubleGridProperty("PORO").getData(), compressedToCartesianElemIdx);

        this->compressed_permx = compressed_copy(ecl3dProps.getDoubleGridProperty("PERMX").getData(), compressedToCartesianElemIdx);

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMY"))
            this->compressed_permy = compressed_copy(ecl3dProps.getDoubleGridProperty("PERMY").getData(), compressedToCartesianElemIdx);
        else
            this->compressed_permy = compressed_copy(ecl3dProps.getDoubleGridProperty("PERMX").getData(), compressedToCartesianElemIdx);

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMZ"))
            this->compressed_permz = compressed_copy(ecl3dProps.getDoubleGridProperty("PERMZ").getData(), compressedToCartesianElemIdx);
        else
            this->compressed_permz = compressed_copy(ecl3dProps.getDoubleGridProperty("PERMZ").getData(), compressedToCartesianElemIdx);
    }

#endif



    unsigned satRegion(std::size_t active_index) const {
        return this->compressed_satnum->operator[](active_index) - 1;
    }

    double permx(std::size_t active_index) const {
        return this->compressed_permx->operator[](active_index);
    }

    double permy(std::size_t active_index) const {
        return this->compressed_permy->operator[](active_index);
    }

    double permz(std::size_t active_index) const {
        return this->compressed_permz->operator[](active_index);
    }

    double poro(std::size_t active_index) const {
        return this->compressed_poro->operator[](active_index);
    }

    const double * swl(std::size_t active_index) const {
        return this->satfunc(this->compressed_swl, active_index);
    }

    const double * sgl(std::size_t active_index) const {
        return this->satfunc(this->compressed_sgl, active_index);
    }

    const double * swcr(std::size_t active_index) const {
        return this->satfunc(this->compressed_swcr, active_index);
    }

    const double * sgcr(std::size_t active_index) const {
        return this->satfunc(this->compressed_sgcr, active_index);
    }

    const double * sowcr(std::size_t active_index) const {
        return this->satfunc(this->compressed_sowcr, active_index);
    }

    const double * sogcr(std::size_t active_index) const {
        return this->satfunc(this->compressed_sogcr, active_index);
    }

    const double * swu(std::size_t active_index) const {
        return this->satfunc(this->compressed_swu, active_index);
    }

    const double * sgu(std::size_t active_index) const {
        return this->satfunc(this->compressed_sgu, active_index);
    }

    const double * pcw(std::size_t active_index) const {
        return this->satfunc(this->compressed_pcw, active_index);
    }

    const double * pcg(std::size_t active_index) const {
        return this->satfunc(this->compressed_pcg, active_index);
    }

    const double * krw(std::size_t active_index) const {
        return this->satfunc(this->compressed_krw, active_index);
    }

    const double * krg(std::size_t active_index) const {
        return this->satfunc(this->compressed_krg, active_index);
    }

    const double * kro(std::size_t active_index) const {
        return this->satfunc(this->compressed_kro, active_index);
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

    const double * satfunc(const std::unique_ptr<std::vector<double>>& data, std::size_t active_index) const {
        if (data)
            return &(data->operator[](active_index));
        return nullptr;
    }


    std::unique_ptr<std::vector<int>> compressed_satnum;
    std::unique_ptr<std::vector<double>> compressed_swl;
    std::unique_ptr<std::vector<double>> compressed_sgl;
    std::unique_ptr<std::vector<double>> compressed_swcr;
    std::unique_ptr<std::vector<double>> compressed_sgcr;
    std::unique_ptr<std::vector<double>> compressed_sowcr;
    std::unique_ptr<std::vector<double>> compressed_sogcr;
    std::unique_ptr<std::vector<double>> compressed_swu;
    std::unique_ptr<std::vector<double>> compressed_sgu;
    std::unique_ptr<std::vector<double>> compressed_pcw;
    std::unique_ptr<std::vector<double>> compressed_pcg;
    std::unique_ptr<std::vector<double>> compressed_krw;
    std::unique_ptr<std::vector<double>> compressed_kro;
    std::unique_ptr<std::vector<double>> compressed_krg;

    std::unique_ptr<std::vector<double>> compressed_permx;
    std::unique_ptr<std::vector<double>> compressed_permy;
    std::unique_ptr<std::vector<double>> compressed_permz;
    std::unique_ptr<std::vector<double>> compressed_poro;
};
}
#endif

