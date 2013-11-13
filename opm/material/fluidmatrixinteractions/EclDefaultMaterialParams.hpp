// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2013 by Andreas Lauser                                    *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Opm::EclDefaultMaterialParams
 */
#ifndef OPM_ECL_DEFAULT_MATERIAL_PARAMS_HPP
#define OPM_ECL_DEFAULT_MATERIAL_PARAMS_HPP

#include <type_traits>

namespace Opm {

/*!
 * \brief Default implementation for the parameters required by the
 *        default three-phase capillary pressure model used by
 *        Eclipse.
 *
 * Essentially, this class just stores the two parameter objects for
 * the twophase capillary pressure laws.
 */
template<class GasOilParams, class OilWaterParams>
class EclDefaultMaterialParams
{
public:
    static_assert(GasOilParams::numPhases == 2,
                  "The number of phases considered by the gas-oil capillary "
                  "pressure law must be two!")
    static_assert(OilWaterParams::numPhases == 2,
                  "The number of phases considered by the oil-water capillary "
                  "pressure law must be two!")
    static_assert(std::is_same<typename GasOilParams::Scalar,
                               typename OilWaterParams::Scalar>::value,
                  "The two two-phase capillary pressure laws must use the same "
                  "type of floating point values.");

    typedef typename GasOilParams::Scalar Scalar;
    enum { numPhases = 3 };

    /*!
     * \brief The default constructor.
     */
    EclDefaultMaterialParams()
    { }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
    const GasOilParams& gasOilParams() const
    { return gasOilParams_; }

    /*!
     * \brief Set the parameter object for the gas-oil twophase law.
     */
    void setGasOilParams(const GasOilParams& val)
    { gasOilParams_ = val; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    const OilWaterParams& oilWaterParams() const
    { return oilWaterParams_; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    void oilWaterParams(const OilWaterParams& val)
    { oilWaterParams_ = val; }

private:
    GasOilParams gasOilParams_;
    OilWaterParams oilWaterParams_;
};
} // namespace Opm

#endif
