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
 * \copydoc Opm::EclTwoPhaseMaterialParams
 */
#ifndef OPM_ECL_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_ECL_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <type_traits>
#include <cassert>
#include <memory>

#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm {
enum class EclTwoPhaseApproach {
    EclTwoPhaseGasOil,
    EclTwoPhaseOilWater,
    EclTwoPhaseGasWater
};

/*!
 * \brief Implementation for the parameters required by the material law for two-phase
 *        simulations.
 *
 * Essentially, this class just stores the two parameter objects for
 * the twophase capillary pressure laws.
 */
template<class Traits, class GasOilParamsT, class OilWaterParamsT, class GasWaterParamsT>
class EclTwoPhaseMaterialParams : public EnsureFinalized
{
    typedef typename Traits::Scalar Scalar;
    enum { numPhases = 3 };
public:
    using EnsureFinalized :: finalize;


    typedef GasOilParamsT GasOilParams;
    typedef OilWaterParamsT OilWaterParams;
    typedef GasWaterParamsT GasWaterParams;

    /*!
     * \brief The default constructor.
     */
    EclTwoPhaseMaterialParams()
    {
    }

    void setApproach(EclTwoPhaseApproach newApproach)
    { approach_ = newApproach; }

    EclTwoPhaseApproach approach() const
    { return approach_; }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
    const GasOilParams& gasOilParams() const
    { EnsureFinalized::check(); return *gasOilParams_; }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
    GasOilParams& gasOilParams()
    { EnsureFinalized::check(); return *gasOilParams_; }

    /*!
     * \brief Set the parameter object for the gas-oil twophase law.
     */
    void setGasOilParams(std::shared_ptr<GasOilParams> val)
    { gasOilParams_ = val; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    const OilWaterParams& oilWaterParams() const
    { EnsureFinalized::check(); return *oilWaterParams_; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    OilWaterParams& oilWaterParams()
    { EnsureFinalized::check(); return *oilWaterParams_; }

    /*!
     * \brief Set the parameter object for the oil-water twophase law.
     */
    void setOilWaterParams(std::shared_ptr<OilWaterParams> val)
    { oilWaterParams_ = val; }

  /*!
     * \brief The parameter object for the gas-water twophase law.
     */
    const GasWaterParams& gasWaterParams() const
    { EnsureFinalized::check(); return *gasWaterParams_; }

    /*!
     * \brief The parameter object for the gas-water twophase law.
     */
    GasWaterParams& gasWaterParams()
    { EnsureFinalized::check(); return *gasWaterParams_; }

    /*!
     * \brief Set the parameter object for the gas-water twophase law.
     */
    void setGasWaterParams(std::shared_ptr<GasWaterParams> val)
    { gasWaterParams_ = val; }
    
private:
    EclTwoPhaseApproach approach_;

    std::shared_ptr<GasOilParams> gasOilParams_;
    std::shared_ptr<OilWaterParams> oilWaterParams_;
    std::shared_ptr<GasWaterParams> gasWaterParams_;
};
} // namespace Opm

#endif
