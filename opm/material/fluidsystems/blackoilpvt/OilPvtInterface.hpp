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
 * \copydoc Opm::OilPvtInterface
 */
#ifndef OPM_OIL_PVT_INTERFACE_HPP
#define OPM_OIL_PVT_INTERFACE_HPP

#include <opm/material/common/OpmFinal.hpp>

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil
 *        phase in the black-oil model.
 *
 * This is the base class which which provides an API for the actual PVT implementation
 * classes which based on dynamic polymorphism. The rationale to use dynamic polymorphism
 * here is that this enables the fluid system to easily switch the used PVT relations for
 * the individual fluid phases.
 *
 * Note that, since the application for this class is the black-oil fluid system, the API
 * exposed by this class is pretty specific to the black-oil model.
 */
template <class Scalar, class Evaluation = Scalar>
class OilPvtInterface
{
public:
    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    virtual Evaluation viscosity(int regionIdx,
                                 const Evaluation& temperature,
                                 const Evaluation& pressure,
                                 const Evaluation& XoG) const = 0;
    virtual Scalar viscosity(int regionIdx,
                             Scalar temperature,
                             Scalar pressure,
                             Scalar XoG) const = 0;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    virtual Evaluation formationVolumeFactor(int regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& XoG) const = 0;
    virtual Scalar formationVolumeFactor(int regionIdx,
                                         Scalar temperature,
                                         Scalar pressure,
                                         Scalar XoG) const = 0;

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    virtual Evaluation density(int regionIdx,
                               const Evaluation& temperature,
                               const Evaluation& pressure,
                               const Evaluation& XoG) const = 0;
    virtual Scalar density(int regionIdx,
                           Scalar temperature,
                           Scalar pressure,
                           Scalar XoG) const = 0;

    /*!
     * \brief Returns the fugacity coefficient [Pa] of a component in the fluid phase given
     *        a set of parameters.
     */
    virtual Evaluation fugacityCoefficient(int regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure,
                                           int compIdx) const = 0;
    virtual Scalar fugacityCoefficient(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure,
                                       int compIdx) const = 0;

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of saturated oil.
     */
    virtual Evaluation gasDissolutionFactor(int regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure) const = 0;
    virtual Scalar gasDissolutionFactor(int regionIdx,
                                        Scalar temperature,
                                        Scalar pressure) const = 0;

    /*!
     * \brief Returns the saturation pressure [Pa] of oil given the mass fraction of the
     *        gas component in the oil phase.
     *
     * Calling this method only makes sense for live oil. All other implementations of
     * the black-oil PVT interface will just throw an exception...
     */
    virtual Evaluation oilSaturationPressure(int regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& XoG) const = 0;
    virtual Scalar oilSaturationPressure(int regionIdx,
                                         Scalar temperature,
                                         Scalar XoG) const = 0;

    /*!
     * \brief Returns the gas mass fraction of gas-saturated oil at a given temperatire
     *        and pressure [-].
     *
     * Calling this method only makes sense for oil. For all other phases an exception
     * will be thrown...
     */
    virtual Evaluation saturatedOilGasMassFraction(int regionIdx,
                                                   const Evaluation& temperature,
                                                   const Evaluation& pressure) const = 0;
    virtual Scalar saturatedOilGasMassFraction(int regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const = 0;

    /*!
     * \brief Returns the gas mole fraction of gas-saturated oil at a given temperatire
     *        and pressure [-].
     *
     * Calling this method only makes sense for oil. For all other phases an exception
     * will be thrown...
     */
    virtual Evaluation saturatedOilGasMoleFraction(int regionIdx,
                                                   const Evaluation& temperature,
                                                   const Evaluation& pressure) const = 0;
    virtual Scalar saturatedOilGasMoleFraction(int regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const = 0;
};

// if Evaluation == Scalar, the PVT interface provides only the scalar variants of the
// methods.
template <class Scalar>
class OilPvtInterface<Scalar, Scalar>
{
public:
    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    virtual Scalar viscosity(int regionIdx,
                             Scalar temperature,
                             Scalar pressure,
                             Scalar XoG) const = 0;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    virtual Scalar formationVolumeFactor(int regionIdx,
                                         Scalar temperature,
                                         Scalar pressure,
                                         Scalar XoG) const = 0;

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    virtual Scalar density(int regionIdx,
                           Scalar temperature,
                           Scalar pressure,
                           Scalar XoG) const = 0;

    /*!
     * \brief Returns the fugacity coefficient [Pa] of a component in the fluid phase given
     *        a set of parameters.
     */
    virtual Scalar fugacityCoefficient(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure,
                                       int compIdx) const = 0;

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of saturated oil.
     */
    virtual Scalar gasDissolutionFactor(int regionIdx,
                                        Scalar temperature,
                                        Scalar pressure) const = 0;

    /*!
     * \brief Returns the saturation pressure [Pa] of oil given the mass fraction of the
     *        gas component in the oil phase.
     *
     * Calling this method only makes sense for live oil. All other implementations of
     * the black-oil PVT interface will just throw an exception...
     */
    virtual Scalar oilSaturationPressure(int regionIdx,
                                         Scalar temperature,
                                         Scalar XoG) const = 0;

    /*!
     * \brief Returns the gas mass fraction of gas-saturated oil at a given temperatire
     *        and pressure [-].
     *
     * Calling this method only makes sense for oil. For all other phases an exception
     * will be thrown...
     */
    virtual Scalar saturatedOilGasMassFraction(int regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const = 0;

    /*!
     * \brief Returns the gas mole fraction of gas-saturated oil at a given temperatire
     *        and pressure [-].
     *
     * Calling this method only makes sense for oil. For all other phases an exception
     * will be thrown...
     */
    virtual Scalar saturatedOilGasMoleFraction(int regionIdx,
                                               Scalar temperature,
                                               Scalar pressure) const = 0;
};

// To prevent the need for most code duplication, derived classes can use this class to
// call templated variants of all the methods
template <class Scalar, class Evaluation, class Implementation>
class OilPvtInterfaceTemplateWrapper
    : public OilPvtInterface<Scalar, Evaluation>
{
    Evaluation viscosity(int regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& XoG) const OPM_FINAL
    { return asImp_().template viscosity_<Evaluation>(regionIdx, temperature, pressure, XoG); }
    Scalar viscosity(int regionIdx,
                     Scalar temperature,
                     Scalar pressure,
                     Scalar XoG) const OPM_FINAL
    { return asImp_().template viscosity_<Scalar>(regionIdx, temperature, pressure, XoG); }

    Evaluation formationVolumeFactor(int regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& pressure,
                                     const Evaluation& XoG) const OPM_FINAL
    { return asImp_().template formationVolumeFactor_<Evaluation>(regionIdx, temperature, pressure, XoG); }
    Scalar formationVolumeFactor(int regionIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 Scalar XoG) const OPM_FINAL
    { return asImp_().template formationVolumeFactor_<Scalar>(regionIdx, temperature, pressure, XoG); }

    Evaluation density(int regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure,
                       const Evaluation& XoG) const OPM_FINAL
    { return asImp_().template density_<Evaluation>(regionIdx, temperature, pressure, XoG); }
    Scalar density(int regionIdx,
                   Scalar temperature,
                   Scalar pressure,
                   Scalar XoG) const OPM_FINAL
    { return asImp_().template density_<Scalar>(regionIdx, temperature, pressure, XoG); }

    Evaluation fugacityCoefficient(int regionIdx,
                                   const Evaluation& temperature,
                                   const Evaluation& pressure,
                                   int compIdx) const OPM_FINAL
    { return asImp_().template fugacityCoefficient_<Evaluation>(regionIdx, temperature, pressure, compIdx); }
    Scalar fugacityCoefficient(int regionIdx,
                               Scalar temperature,
                               Scalar pressure,
                               int compIdx) const OPM_FINAL
    { return asImp_().template fugacityCoefficient_<Scalar>(regionIdx, temperature, pressure, compIdx); }

    Evaluation gasDissolutionFactor(int regionIdx,
                                    const Evaluation& temperature,
                                    const Evaluation& pressure) const OPM_FINAL
    { return asImp_().template gasDissolutionFactor_<Evaluation>(regionIdx, temperature, pressure); }
    Scalar gasDissolutionFactor(int regionIdx,
                                Scalar temperature,
                                Scalar pressure) const OPM_FINAL
    { return asImp_().template gasDissolutionFactor_<Scalar>(regionIdx, temperature, pressure); }

    Evaluation oilSaturationPressure(int regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& XoG) const OPM_FINAL
    { return asImp_().template oilSaturationPressure_<Evaluation>(regionIdx, temperature, XoG); }
    Scalar oilSaturationPressure(int regionIdx,
                                 Scalar temperature,
                                 Scalar XoG) const OPM_FINAL
    { return asImp_().template oilSaturationPressure_<Scalar>(regionIdx, temperature, XoG); }

    Evaluation saturatedOilGasMassFraction(int regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure) const OPM_FINAL
    { return asImp_().template saturatedOilGasMassFraction_<Evaluation>(regionIdx, temperature, pressure); }
    Scalar saturatedOilGasMassFraction(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure) const OPM_FINAL
    { return asImp_().template saturatedOilGasMassFraction_<Scalar>(regionIdx, temperature, pressure); }

    Evaluation saturatedOilGasMoleFraction(int regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure) const OPM_FINAL
    { return asImp_().template saturatedOilGasMoleFraction_<Evaluation>(regionIdx, temperature, pressure); }
    Scalar saturatedOilGasMoleFraction(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure) const OPM_FINAL
    { return asImp_().template saturatedOilGasMoleFraction_<Scalar>(regionIdx, temperature, pressure); }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

template <class Scalar, class Implementation>
class OilPvtInterfaceTemplateWrapper<Scalar, Scalar, Implementation>
    : public OilPvtInterface<Scalar, Scalar>
{
    Scalar viscosity(int regionIdx,
                     Scalar temperature,
                     Scalar pressure,
                     Scalar XoG) const OPM_FINAL
    { return asImp_().template viscosity_<Scalar>(regionIdx, temperature, pressure, XoG); }
    Scalar formationVolumeFactor(int regionIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 Scalar XoG) const OPM_FINAL
    { return asImp_().template formationVolumeFactor_<Scalar>(regionIdx, temperature, pressure, XoG); }
    Scalar density(int regionIdx,
                   Scalar temperature,
                   Scalar pressure,
                   Scalar XoG) const OPM_FINAL
    { return asImp_().template density_<Scalar>(regionIdx, temperature, pressure, XoG); }
    Scalar fugacityCoefficient(int regionIdx,
                               Scalar temperature,
                               Scalar pressure,
                               int compIdx) const OPM_FINAL
    { return asImp_().template fugacityCoefficient_<Scalar>(regionIdx, temperature, pressure, compIdx); }
    Scalar gasDissolutionFactor(int regionIdx,
                                Scalar temperature,
                                Scalar pressure) const OPM_FINAL
    { return asImp_().template gasDissolutionFactor_<Scalar>(regionIdx, temperature, pressure); }
    Scalar oilSaturationPressure(int regionIdx,
                                 Scalar temperature,
                                 Scalar XoG) const OPM_FINAL
    { return asImp_().template oilSaturationPressure_<Scalar>(regionIdx, temperature, XoG); }
    Scalar saturatedOilGasMassFraction(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure) const OPM_FINAL
    { return asImp_().template saturatedOilGasMassFraction_<Scalar>(regionIdx, temperature, pressure); }
    Scalar saturatedOilGasMoleFraction(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure) const OPM_FINAL
    { return asImp_().template saturatedOilGasMoleFraction_<Scalar>(regionIdx, temperature, pressure); }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
};


} // namespace Opm

#endif
