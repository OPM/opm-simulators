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
 * \copydoc Opm::EclThconrLawParams
 */
#ifndef OPM_ECL_THCONR_LAW_PARAMS_HPP
#define OPM_ECL_THCONR_LAW_PARAMS_HPP

#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm {

/*!
 * \brief The default implementation of a parameter object for the
 *        thermal conduction law based on the THCONR keyword from ECL.
 */
template <class ScalarT>
class EclThconrLawParams : public EnsureFinalized
{
public:
    typedef ScalarT Scalar;

    EclThconrLawParams(const EclThconrLawParams&) = default;

    EclThconrLawParams()
    { }

    /*!
     * \brief Set the total thermal conductivity [J/m^2 / (K/m)] of at Sg = 0
     */
    void setReferenceTotalThermalConductivity(Scalar value)
    { referenceTotalThermalConductivity_ = value; }

    /*!
     * \brief The total thermal conductivity [J/m^2 / (K/m)] of at Sg = 0
     */
    Scalar referenceTotalThermalConductivity() const
    { EnsureFinalized::check(); return referenceTotalThermalConductivity_; }

    /*!
     * \brief Set the gas saturation dependence of thermal conductivity [-]
     */
    void setDTotalThermalConductivity_dSg(Scalar value)
    { dTotalThermalConductivity_dSg_ = value; }

    /*!
     * \brief The gas saturation dependence of thermal conductivity [-]
     */
    Scalar dTotalThermalConductivity_dSg() const
    { EnsureFinalized::check(); return dTotalThermalConductivity_dSg_; }

private:
    Scalar referenceTotalThermalConductivity_;
    Scalar dTotalThermalConductivity_dSg_;
};

} // namespace Opm

#endif
