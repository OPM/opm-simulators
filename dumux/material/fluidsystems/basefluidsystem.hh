// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Dumux::BaseFluidSystem
 */
#ifndef DUMUX_BASE_FLUID_SYSTEM_HH
#define DUMUX_BASE_FLUID_SYSTEM_HH

#include "nullparametercache.hh"

#include <dumux/common/exceptions.hh>

#include <dune/common/classname.hh>

namespace Dumux {

/*!
 * \ingroup Fluidsystems
 * \brief The base class for all fluid systems.
 */
template <class Scalar, class Implementation>
class BaseFluidSystem
{
public:
    /*!
     * \brief The type of the fluid system's parameter cache
     *
     * The parameter cache can be used to avoid re-calculating
     * expensive parameters for multiple quantities. Be aware that
     * what the parameter cache actually does is specific for each
     * fluid system and that it is opaque outside the fluid system.
     */
    typedef NullParameterCache ParameterCache;

    //! Number of chemical species in the fluid system
    static const int numComponents = -1000;

    //! Number of fluid phases in the fluid system
    static const int numPhases = -2000;

    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    static char *phaseName(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The fluid system '" << Dune::className<Implementation>() << "' does not provide a phaseName() method!");
    }

    /*!
     * \brief Return whether a phase is liquid
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    static bool isLiquid(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The fluid system '" << Dune::className<Implementation>() << "' does not provide a isLiquid() method!");
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if Henry's law and Rault's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    static bool isIdealMixture(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The fluid system '" << Dune::className<Implementation>() << "' does not provide a isIdealMixture() method!");
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    static bool isCompressible(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The fluid system '" << Dune::className<Implementation>() << "' does not provide a isCompressible() method!");
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    static bool isIdealGas(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The fluid system '" << Dune::className<Implementation>() << "' does not provide a isIdealGas() method!");
    }

    /*!
     * \brief Return the human readable name of a component
     *
     * \copydoc Doxygen::compIdxParam
     */
    static const char *componentName(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The fluid system '" << Dune::className<Implementation>() << "' does not provide a componentName() method!");
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     *
     * \copydoc Doxygen::compIdxParam
     */
    static Scalar molarMass(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The fluid system '" << Dune::className<Implementation>() << "' does not provide a molarMass() method!");
    }

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    { }

    /*!
     * \brief Calculate the density [kg/m^3] of a fluid phase
     *
     * \copydoc Doxygen::fluidSystemBaseParams
     * \copydoc Doxygen::phaseIdxParam
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
     * \f[ f_\kappa = \phi_\kappa\,x_{\kappa} \f]
     *
     * \copydoc Doxygen::fluidSystemBaseParams
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
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
     *
     * \copydoc Doxygen::fluidSystemBaseParams
     * \copydoc Doxygen::phaseIdxParam
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
     * gradient of the mole fraction and follows the law
     *
     * \f[ J = - D \mathbf{grad} x^\kappa_\alpha \f]
     *
     * where \f$x_\alpha^\kappa\f$ is the component's mole fraction in
     * phase \f$\alpha\f$, \f$D\f$ is the diffusion coefficient and
     * \f$J\f$ is the diffusive flux.
     *
     * \copydoc Doxygen::fluidSystemBaseParams
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
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
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy [J/kg].
     *
     * \copydoc Doxygen::fluidSystemBaseParams
     * \copydoc Doxygen::phaseIdxParam
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
     * \copydoc Doxygen::fluidSystemBaseParams
     * \copydoc Doxygen::phaseIdxParam
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
     * \copydoc Doxygen::fluidSystemBaseParams
     * \copydoc Doxygen::phaseIdxParam
     */
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
