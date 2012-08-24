// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Christoph Grueninger                              *
 *   Copyright (C) 2012 by Philipp Nuske                                     *
 *   Copyright (C) 2012 by Vishal Jambhekar                                  *
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
 * \ingroup Fluidsystems
 *
 * \brief Fluid system base class.
 */
#ifndef DUMUX_BASE_FLUID_SYSTEM_HH
#define DUMUX_BASE_FLUID_SYSTEM_HH

#include <dune/common/classname.hh>

#include <dumux/common/exceptions.hh>

namespace Dumux
{
/*!
 * \ingroup Fluidsystems
 * \brief Fluid system base class.
 */
template <class Scalar, class Implementation>
class BaseFluidSystem
{
public:
    /*!
     * \brief Calculate the density [kg/m^3] of a fluid phase
     */
    template <class FluidState, class ParameterCache>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The fluid system '" << Dune::className<Implementation>() << "' does not provide a density() method!");
    }

    /*!
     * \brief Calculate the fugacity coefficient [Pa] of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\phi_\kappa\f$ is connected to the
     * fugacity \f$f_\kappa\f$ and the component's molarity
     * \f$x_\kappa\f$ by means of the relation
     *
     * \f[ f_\kappa = \phi_\kappa * x_{\kappa} \f]
     */
    template <class FluidState, class ParameterCache>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "The fluid system '" << Dune::className<Implementation>() << "'  does not provide a fugacityCoefficient() method!");
    }

    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s]
     */
    template <class FluidState, class ParameterCache>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "The fluid system '" << Dune::className<Implementation>() << "'  does not provide a viscosity() method!");
    }

    /*!
     * \brief Calculate the binary molecular diffusion coefficient for
     *        a component in a fluid phase [mol^2 * s / (kg*m^3)]
     *
     * Molecular diffusion of a compoent \f$\kappa\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \mathbf{grad} \mu_\kappa \f]
     *
     * where \f$\mu_\kappa\f$ is the component's chemical potential,
     * \f$D\f$ is the diffusion coefficient and \f$J\f$ is the
     * diffusive flux. \f$\mu_\kappa\f$ is connected to the component's
     * fugacity \f$f_\kappa\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$p_\alpha\f$ and \f$T_\alpha\f$ are the fluid phase'
     * pressure and temperature.
     */
    template <class FluidState, class ParameterCache>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    {
         DUNE_THROW(Dune::NotImplemented, "The fluid system '" << Dune::className<Implementation>() << "'  does not provide a diffusionCoefficient() method!");
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     */
    template <class FluidState, class ParameterCache>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             const ParameterCache &paramCache,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        DUNE_THROW(Dune::NotImplemented, "The fluid system '" << Dune::className<Implementation>() << "'  does not provide a binaryDiffusionCoefficient() method!");
    }

    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy [J/kg].
     *
     *  \todo This fluid system neglects the contribution of
     *        gas-molecules in the liquid phase. This contribution is
     *        probably not big. Somebody would have to find out the
     *        enthalpy of solution for this system. ...
     */
    template <class FluidState, class ParameterCache>
    static Scalar enthalpy(const FluidState &fluidState,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "The fluid system '" << Dune::className<Implementation>() << "'  does not provide an enthalpy() method!");
    }

    /*!
     * \brief Thermal conductivity of a fluid phase [W/(m K)].
     *
     * Use the conductivity of air and water as a first approximation.
     * Source:
     * http://en.wikipedia.org/wiki/List_of_thermal_conductivities
     */
    template <class FluidState, class ParameterCache>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "The fluid system '" << Dune::className<Implementation>() << "'  does not provide a thermalConductivity() method!");
    }

    /*!
     * \brief Specific isobaric heat capacity of a fluid phase [J/kg].
     *
     * \param paramCache   mutable parameters
     * \param phaseIdx  for which phase to give back the heat capacity
     * \param fluidState represents all relevant thermodynamic quantities of a
     *  fluid system
     * */
    template <class FluidState, class ParameterCache>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "The fluid system '" << Dune::className<Implementation>() << "'  does not provide a heatCapacity() method!");
    }
};

} // end namepace

#endif
