// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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
 * \copydoc Opm::GasPvtInterface
 */
#ifndef OPM_GAS_PVT_INTERFACE_HPP
#define OPM_GAS_PVT_INTERFACE_HPP

#include <opm/material/common/OpmFinal.hpp>

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas
 *        phase in the black-oil model.
 *
 * This is the base class which which provides an API for the actual PVT implementation
 * classes which based on dynamic polymorphism. The rationale to use dynamic polymorphism
 * here is that this enables the fluid system to easily switch the used PVT relations for
 * the individual fluid phases.
 *
 * Note that, since the application for this class is the black oil fluid system, the API
 * exposed by this class is pretty specific to the black oil model.
 */
template <class Scalar, class Evaluation = Scalar>
class GasPvtInterface
{
public:
    virtual ~GasPvtInterface() {}

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    virtual Evaluation viscosity(unsigned regionIdx,
                                 const Evaluation& temperature,
                                 const Evaluation& pressure,
                                 const Evaluation& XgO) const = 0;
    virtual Scalar viscosity(unsigned regionIdx,
                             Scalar temperature,
                             Scalar pressure,
                             Scalar XgO) const = 0;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    virtual Evaluation formationVolumeFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& XgO) const = 0;
    virtual Scalar formationVolumeFactor(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar pressure,
                                         Scalar XgO) const = 0;

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    virtual Evaluation density(unsigned regionIdx,
                               const Evaluation& temperature,
                               const Evaluation& pressure,
                               const Evaluation& XgO) const = 0;
    virtual Scalar density(unsigned regionIdx,
                           Scalar temperature,
                           Scalar pressure,
                           Scalar XgO) const = 0;

    /*!
     * \brief Returns the fugacity coefficient [Pa] of a component in the fluid phase given
     *        a set of parameters.
     */
    virtual Evaluation fugacityCoefficient(unsigned regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure,
                                           unsigned compIdx) const = 0;
    virtual Scalar fugacityCoefficient(unsigned regionIdx,
                                       Scalar temperature,
                                       Scalar pressure,
                                       unsigned compIdx) const = 0;

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of oil saturated gas.
     */
    virtual Evaluation oilVaporizationFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const = 0;
    virtual Scalar oilVaporizationFactor(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar pressure) const = 0;

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param XgO The mass fraction of the oil component in the gas phase [-]
     */
    virtual Evaluation gasSaturationPressure(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& XgO) const = 0;
    virtual Scalar gasSaturationPressure(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar XgO) const = 0;

    /*!
     * \brief Returns the gas mass fraction of oil-saturated gas at a given temperatire
     *        and pressure [-].
     */
    virtual Evaluation saturatedGasOilMassFraction(unsigned regionIdx,
                                                   const Evaluation& temperature,
                                                   const Evaluation& pressure) const = 0;
    virtual Scalar saturatedGasOilMassFraction(unsigned regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const = 0;

    /*!
     * \brief Returns the gas mole fraction of oil-saturated gas at a given temperatire
     *        and pressure [-].
     */
    virtual Evaluation saturatedGasOilMoleFraction(unsigned regionIdx,
                                                   const Evaluation& temperature,
                                                   const Evaluation& pressure) const = 0;
    virtual Scalar saturatedGasOilMoleFraction(unsigned regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const = 0;
};

template <class Scalar>
class GasPvtInterface<Scalar, Scalar>
{
public:
    virtual Scalar viscosity(unsigned regionIdx,
                             Scalar temperature,
                             Scalar pressure,
                             Scalar XgO) const = 0;
    virtual Scalar formationVolumeFactor(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar pressure,
                                         Scalar XgO) const = 0;
    virtual Scalar density(unsigned regionIdx,
                           Scalar temperature,
                           Scalar pressure,
                           Scalar XgO) const = 0;
    virtual Scalar fugacityCoefficient(unsigned regionIdx,
                                       Scalar temperature,
                                       Scalar pressure,
                                       unsigned compIdx) const = 0;
    virtual Scalar oilVaporizationFactor(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar pressure) const = 0;
    virtual Scalar gasSaturationPressure(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar XgO) const = 0;
    virtual Scalar saturatedGasOilMassFraction(unsigned regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const = 0;
    virtual Scalar saturatedGasOilMoleFraction(unsigned regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const = 0;
};

template <class Scalar, class Evaluation, class Implementation>
class GasPvtInterfaceTemplateWrapper : public GasPvtInterface<Scalar, Evaluation>
{
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

public:
    virtual Evaluation viscosity(unsigned regionIdx,
                                 const Evaluation& temperature,
                                 const Evaluation& pressure,
                                 const Evaluation& XgO) const OPM_FINAL
    { return asImp_().viscosity_(regionIdx, temperature, pressure, XgO); }
    virtual Scalar viscosity(unsigned regionIdx,
                             Scalar temperature,
                             Scalar pressure,
                             Scalar XgO) const OPM_FINAL
    { return asImp_().viscosity_(regionIdx, temperature, pressure, XgO); }

    virtual Evaluation formationVolumeFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& XgO) const OPM_FINAL
    { return asImp_().formationVolumeFactor_(regionIdx, temperature, pressure, XgO); }
    virtual Scalar formationVolumeFactor(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar pressure,
                                         Scalar XgO) const OPM_FINAL
    { return asImp_().formationVolumeFactor_(regionIdx, temperature, pressure, XgO); }

    virtual Evaluation density(unsigned regionIdx,
                               const Evaluation& temperature,
                               const Evaluation& pressure,
                               const Evaluation& XgO) const OPM_FINAL
    { return asImp_().density_(regionIdx, temperature, pressure, XgO); }
    virtual Scalar density(unsigned regionIdx,
                           Scalar temperature,
                           Scalar pressure,
                           Scalar XgO) const OPM_FINAL
    { return asImp_().density_(regionIdx, temperature, pressure, XgO); }

    virtual Evaluation fugacityCoefficient(unsigned regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure,
                                           unsigned compIdx) const OPM_FINAL
    { return asImp_().fugacityCoefficient_(regionIdx, temperature, pressure, compIdx); }
    virtual Scalar fugacityCoefficient(unsigned regionIdx,
                                       Scalar temperature,
                                       Scalar pressure,
                                       unsigned compIdx) const OPM_FINAL
    { return asImp_().fugacityCoefficient_(regionIdx, temperature, pressure, compIdx); }

    virtual Evaluation oilVaporizationFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const OPM_FINAL
    { return asImp_().oilVaporizationFactor_(regionIdx, temperature, pressure); }
    virtual Scalar oilVaporizationFactor(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar pressure) const OPM_FINAL
    { return asImp_().oilVaporizationFactor_(regionIdx, temperature, pressure); }

    virtual Evaluation gasSaturationPressure(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& XgO) const OPM_FINAL
    { return asImp_().gasSaturationPressure_(regionIdx, temperature, XgO); }
    virtual Scalar gasSaturationPressure(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar XgO) const OPM_FINAL
    { return asImp_().gasSaturationPressure_(regionIdx, temperature, XgO); }

    virtual Evaluation saturatedGasOilMassFraction(unsigned regionIdx,
                                                   const Evaluation& temperature,
                                                   const Evaluation& pressure) const OPM_FINAL
    { return asImp_().saturatedGasOilMassFraction_(regionIdx, temperature, pressure); }
    virtual Scalar saturatedGasOilMassFraction(unsigned regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const OPM_FINAL
    { return asImp_().saturatedGasOilMassFraction_(regionIdx, temperature, pressure); }

    virtual Evaluation saturatedGasOilMoleFraction(unsigned regionIdx,
                                                   const Evaluation& temperature,
                                                   const Evaluation& pressure) const OPM_FINAL
    { return asImp_().saturatedGasOilMoleFraction_(regionIdx, temperature, pressure); }
    virtual Scalar saturatedGasOilMoleFraction(unsigned regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const OPM_FINAL
    { return asImp_().saturatedGasOilMoleFraction_(regionIdx, temperature, pressure); }
};

template <class Scalar, class Implementation>
class GasPvtInterfaceTemplateWrapper<Scalar, Scalar, Implementation>
    : public GasPvtInterface<Scalar, Scalar>
{
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

public:
    virtual Scalar viscosity(unsigned regionIdx,
                             Scalar temperature,
                             Scalar pressure,
                             Scalar XgO) const OPM_FINAL
    { return asImp_().viscosity_(regionIdx, temperature, pressure, XgO); }
    virtual Scalar formationVolumeFactor(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar pressure,
                                         Scalar XgO) const OPM_FINAL
    { return asImp_().formationVolumeFactor_(regionIdx, temperature, pressure, XgO); }
    virtual Scalar density(unsigned regionIdx,
                           Scalar temperature,
                           Scalar pressure,
                           Scalar XgO) const OPM_FINAL
    { return asImp_().density_(regionIdx, temperature, pressure, XgO); }
    virtual Scalar fugacityCoefficient(unsigned regionIdx,
                                       Scalar temperature,
                                       Scalar pressure,
                                       unsigned compIdx) const OPM_FINAL
    { return asImp_().fugacityCoefficient_(regionIdx, temperature, pressure, compIdx); }
    virtual Scalar oilVaporizationFactor(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar pressure) const OPM_FINAL
    { return asImp_().oilVaporizationFactor_(regionIdx, temperature, pressure); }
    virtual Scalar gasSaturationPressure(unsigned regionIdx,
                                         Scalar temperature,
                                         Scalar XgO) const OPM_FINAL
    { return asImp_().gasSaturationPressure_(regionIdx, temperature, XgO); }
    virtual Scalar saturatedGasOilMassFraction(unsigned regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const OPM_FINAL
    { return asImp_().saturatedGasOilMassFraction_(regionIdx, temperature, pressure); }
    virtual Scalar saturatedGasOilMoleFraction(unsigned regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const OPM_FINAL
    { return asImp_().saturatedGasOilMoleFraction_(regionIdx, temperature, pressure); }
};

} // namespace Opm

#endif
