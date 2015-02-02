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
template <class Scalar>
class OilPvtInterface
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

} // namespace Opm

#endif
