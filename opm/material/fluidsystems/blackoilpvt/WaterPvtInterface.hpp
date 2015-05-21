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
 * \copydoc Opm::WaterPvtInterface
 */
#ifndef OPM_WATER_PVT_INTERFACE_HPP
#define OPM_WATER_PVT_INTERFACE_HPP

#include <opm/material/common/OpmFinal.hpp>

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the water
 *        phase in the black-oil model.
 *
 * This is the base class which which provides an API for the actual PVT implementation
 * classes which based on dynamic polymorphism. The rationale to use dynamic polymorphism
 * here is that this enables the fluid system to easily switch the used PVT relations for
 * the individual fluid phases.
 *
 * Note that, since the application for this class is the black oil fluid system, the API
 * exposed by this class is pretty specific to the black-oil model.
 */
template <class Scalar, class Evaluation = Scalar>
class WaterPvtInterface
{
public:
    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    virtual Evaluation viscosity(int regionIdx,
                                 const Evaluation& temperature,
                                 const Evaluation& pressure) const = 0;
    virtual Scalar viscosity(int regionIdx,
                             Scalar temperature,
                             Scalar pressure) const = 0;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    virtual Evaluation formationVolumeFactor(int regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const = 0;
    virtual Scalar formationVolumeFactor(int regionIdx,
                                         Scalar temperature,
                                         Scalar pressure) const = 0;

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    virtual Evaluation density(int regionIdx,
                               const Evaluation& temperature,
                               const Evaluation& pressure) const = 0;
    virtual Scalar density(int regionIdx,
                           Scalar temperature,
                           Scalar pressure) const = 0;

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
};

template <class Scalar>
class WaterPvtInterface<Scalar, Scalar>
{
public:
    virtual Scalar viscosity(int regionIdx,
                             Scalar temperature,
                             Scalar pressure) const = 0;

    virtual Scalar formationVolumeFactor(int regionIdx,
                                         Scalar temperature,
                                         Scalar pressure) const = 0;
    virtual Scalar density(int regionIdx,
                           Scalar temperature,
                           Scalar pressure) const = 0;
    virtual Scalar fugacityCoefficient(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure,
                                       int compIdx) const = 0;
};

template <class Scalar, class Evaluation, class Implementation>
class WaterPvtInterfaceTemplateWrapper : public WaterPvtInterface<Scalar, Evaluation>
{
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

public:
    Evaluation viscosity(int regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure) const OPM_FINAL
    { return asImp_().viscosity_(regionIdx, temperature, pressure); }
    Scalar viscosity(int regionIdx,
                     Scalar temperature,
                     Scalar pressure) const OPM_FINAL
    { return asImp_().viscosity_(regionIdx, temperature, pressure); }

    Evaluation formationVolumeFactor(int regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& pressure) const OPM_FINAL
    { return asImp_().formationVolumeFactor_(regionIdx, temperature, pressure); }
    Scalar formationVolumeFactor(int regionIdx,
                                 Scalar temperature,
                                 Scalar pressure) const OPM_FINAL
    { return asImp_().formationVolumeFactor_(regionIdx, temperature, pressure); }

    Evaluation density(int regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure) const OPM_FINAL
    { return asImp_().density_(regionIdx, temperature, pressure); }
    Scalar density(int regionIdx,
                   Scalar temperature,
                   Scalar pressure) const OPM_FINAL
    { return asImp_().density_(regionIdx, temperature, pressure); }

    Evaluation fugacityCoefficient(int regionIdx,
                                   const Evaluation& temperature,
                                   const Evaluation& pressure,
                                   int compIdx) const OPM_FINAL
    { return asImp_().fugacityCoefficient_(regionIdx, temperature, pressure, compIdx); }
    Scalar fugacityCoefficient(int regionIdx,
                               Scalar temperature,
                               Scalar pressure,
                               int compIdx) const OPM_FINAL
    { return asImp_().fugacityCoefficient_(regionIdx, temperature, pressure, compIdx); }
};

template <class Scalar, class Implementation>
class WaterPvtInterfaceTemplateWrapper<Scalar, Scalar, Implementation>
    : public WaterPvtInterface<Scalar, Scalar>
{
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

public:
    Scalar viscosity(int regionIdx,
                     Scalar temperature,
                     Scalar pressure) const OPM_FINAL
    { return asImp_().viscosity_(regionIdx, temperature, pressure); }

    Scalar formationVolumeFactor(int regionIdx,
                                 Scalar temperature,
                                 Scalar pressure) const OPM_FINAL
    { return asImp_().formationVolumeFactor_(regionIdx, temperature, pressure); }

    Scalar density(int regionIdx,
                   Scalar temperature,
                   Scalar pressure) const OPM_FINAL
    { return asImp_().density_(regionIdx, temperature, pressure); }

    Scalar fugacityCoefficient(int regionIdx,
                               Scalar temperature,
                               Scalar pressure,
                               int compIdx) const OPM_FINAL
    { return asImp_().fugacityCoefficient_(regionIdx, temperature, pressure, compIdx); }
};

} // namespace Opm

#endif
