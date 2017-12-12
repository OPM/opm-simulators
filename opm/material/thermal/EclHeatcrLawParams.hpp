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
 * \copydoc Opm::EclHeatcrLawParams
 */
#ifndef OPM_ECL_HEATCR_LAW_PARAMS_HPP
#define OPM_ECL_HEATCR_LAW_PARAMS_HPP

#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm {

/*!
 * \brief The default implementation of a parameter object for the
 *        ECL thermal law.
 */
template <class ScalarT>
class EclHeatcrLawParams : public EnsureFinalized
{
public:
    typedef ScalarT Scalar;

    EclHeatcrLawParams(const EclHeatcrLawParams&) = default;

    EclHeatcrLawParams()
    { }

    /*!
     * \brief Set the reference temperature for the thermal law.
     *
     * This is a bit hacky because only one temperature is possible, but some
     * memory is saved this way. TODO: Solve this in a better way.
     */
    static void setReferenceTemperature(Scalar value)
    { referenceTemperature_ = value; }

    /*!
     * \brief Return the reference temperature for the thermal law.
     */
    static Scalar referenceTemperature()
    { return referenceTemperature_; }

    /*!
     * \brief Set the specific heat capacity of rock.
     */
    void setReferenceRockHeatCapacity(Scalar value)
    { referenceRockHeatCapacity_ = value; }

    /*!
     * \brief The specific heat capacity of rock.
     */
    Scalar referenceRockHeatCapacity() const
    { EnsureFinalized::check(); return referenceRockHeatCapacity_; }

    /*!
     * \brief Set the derivative of the specific heat capacity of rock w.r.t. temperature.
     */
    void setDRockHeatCapacity_dT(Scalar value)
    { dRockHeatCapacity_dT_ = value; }

    /*!
     * \brief The derivative of the specific heat capacity of rock w.r.t. temperature.
     */
    Scalar dRockHeatCapacity_dT() const
    { EnsureFinalized::check(); return dRockHeatCapacity_dT_; }

private:
    static Scalar referenceTemperature_;

    Scalar referenceRockHeatCapacity_;
    Scalar dRockHeatCapacity_dT_;
};

template <class ScalarT>
ScalarT EclHeatcrLawParams<ScalarT>::referenceTemperature_ = 273.15 + 15.56; // [K]

} // namespace Opm

#endif
