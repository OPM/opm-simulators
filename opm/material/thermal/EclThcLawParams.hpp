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
 * \copydoc Opm::EclThcLawParams
 */
#ifndef OPM_ECL_THC_LAW_PARAMS_HPP
#define OPM_ECL_THC_LAW_PARAMS_HPP

#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm {

/*!
 * \brief The default implementation of a parameter object for the
 *        thermal conduction law based on the THC* keywords from ECL.
 */
template <class ScalarT>
class EclThcLawParams : public EnsureFinalized
{
public:
    typedef ScalarT Scalar;

    EclThcLawParams(const EclThcLawParams&) = default;

    EclThcLawParams()
    { }

    /*!
     * \brief Set the porosity
     */
    void setPorosity(Scalar value)
    { porosity_ = value; }

    /*!
     * \brief Return the porosity
     */
    Scalar porosity() const
    { EnsureFinalized::check(); return porosity_; }

    /*!
     * \brief Set thermal conductivity of pure rock [W/(m*K)]
     */
    void setThcrock(Scalar value)
    { thcrock_ = value; }

    /*!
     * \brief Return thermal conductivity of pure rock [W/(m*K)]
     */
    Scalar thcrock() const
    { EnsureFinalized::check(); return thcrock_; }

    /*!
     * \brief Set thermal conductivity of pure oil [W/(m*K)]
     */
    void setThcoil(Scalar value)
    { thcoil_ = value; }

    /*!
     * \brief Return thermal conductivity of pure oil [W/(m*K)]
     */
    Scalar thcoil() const
    { EnsureFinalized::check(); return thcoil_; }

    /*!
     * \brief Set thermal conductivity of pure gas [W/(m*K)]
     */
    void setThcgas(Scalar value)
    { thcgas_ = value; }

    /*!
     * \brief Return thermal conductivity of pure gas [W/(m*K)]
     */
    Scalar thcgas() const
    { EnsureFinalized::check(); return thcgas_; }

    /*!
     * \brief Set thermal conductivity of pure water [W/(m*K)]
     */
    void setThcwater(Scalar value)
    { thcwater_ = value; }

    /*!
     * \brief Return thermal conductivity of pure water [W/(m*K)]
     */
    Scalar thcwater() const
    { EnsureFinalized::check(); return thcwater_; }

private:
    Scalar porosity_;
    Scalar thcrock_;
    Scalar thcoil_;
    Scalar thcgas_;
    Scalar thcwater_;
};

} // namespace Opm

#endif
