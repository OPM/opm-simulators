/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \brief Specifies the basic API for mutable parameter objects.
 *
 * Fluid systems which do not need any complicated mixture parameters
 * just need this class to be happy.
 */
#ifndef BASIC_MUTABLE_PARAMETERS_HH
#define BASIC_MUTABLE_PARAMETERS_HH

#include <dumux/material/fluidstates/genericfluidstate.hh>
#include <assert.h>

namespace Dumux
{
/*!
 * \brief Specifies the basic API for mutable parameter objects.
 *
 * Fluid systems which do not need any complicated mixture parameters
 * just need this class to be happy.
 */
template <class Scalar,
          class StaticParams,
          class FluidStateT = GenericFluidState<Scalar, StaticParams> >
class BasicMutableParameters
{
    enum { numPhases = StaticParams::numPhases };
    enum { numComponents = StaticParams::numComponents };

public:
    /*!
     * \brief The type of fluid state objects.
     */
    typedef FluidStateT FluidState;

    /*!
     * \brief Returns the fluid state for which the parameter object
     *        is valid.
     */
    FluidState &fluidState() 
    { return fluidState_; }

    /*!
     * \brief Returns the fluid state for which the parameter object
     *        is valid.
     */
    const FluidState &fluidState()  const
    { return fluidState_; }

    /*!
     * \brief Update all parameters of a phase which only depend on 
     *        temperature and/or pressure.
     *
     * This usually means the parameters for the pure components.
     */
    void updatePure(int phaseIdx)
    { }

    /*!
     * \brief Update all parameters of a phase which depend on the
     *        fluid composition. It is assumed that updatePure() has
     *        been called before this method.
     *
     * Here, the mixing rule kicks in.
     */
    void updateMix(int phaseIdx)
    {
        fluidState_.updateMeanMolarMass(phaseIdx);
    };
  
private:
    FluidState fluidState_;
};

} // end namespace Dumux

#endif
