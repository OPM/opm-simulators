/*
  Copyright (C) 2013 by Andreas Lauser

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
 * \copydoc Opm::EclDefaultMaterialParams
 */
#ifndef OPM_ECL_DEFAULT_MATERIAL_PARAMS_HPP
#define OPM_ECL_DEFAULT_MATERIAL_PARAMS_HPP

#include <type_traits>

#include <cassert>

namespace Opm {

/*!
 * \brief Default implementation for the parameters required by the
 *        default three-phase capillary pressure model used by
 *        Eclipse.
 *
 * Essentially, this class just stores the two parameter objects for
 * the twophase capillary pressure laws.
 */
template<class Traits, class GasOilParams, class OilWaterParams>
class EclDefaultMaterialParams
{
    typedef typename Traits::Scalar Scalar;
    enum { numPhases = 3 };
public:
    /*!
     * \brief The default constructor.
     */
    EclDefaultMaterialParams()
    {
#ifndef NDEBUG
        finalized_ = false;
#endif
    }

    /*!
     * \brief Finish the initialization of the parameter object.
     */
    void finalize()
    {
        // Do nothing: The two two-phase parameter objects need to be finalized themselfs!
#ifndef NDEBUG
        finalized_ = true;
#endif
    }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
    const GasOilParams& gasOilParams() const
    { assertFinalized_(); return gasOilParams_; }

    /*!
     * \brief Set the parameter object for the gas-oil twophase law.
     */
    void setGasOilParams(const GasOilParams& val)
    { gasOilParams_ = val; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    const OilWaterParams& oilWaterParams() const
    { assertFinalized_(); return oilWaterParams_; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    void oilWaterParams(const OilWaterParams& val)
    { oilWaterParams_ = val; }

private:
#ifndef NDEBUG
    void assertFinalized_() const
    { assert(finalized_); }

    bool finalized_;
#else
    void assertFinalized_() const
    { }
#endif

    GasOilParams gasOilParams_;
    OilWaterParams oilWaterParams_;
};
} // namespace Opm

#endif
