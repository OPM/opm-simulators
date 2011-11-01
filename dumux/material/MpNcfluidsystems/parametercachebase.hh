/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief The base class of the parameter cache classes for fluid systems
 */
#ifndef DUMUX_PARAMETER_CACHE_BASE_HH
#define DUMUX_PARAMETER_CACHE_BASE_HH

namespace Dumux
{
/*!
 * \brief The base class of the parameter cache classes for fluid systems
 */
template <class Implementation>
class ParameterCacheBase
{
public:
    ParameterCacheBase()
    {};

    template <class FluidState>
    void updateAll(const FluidState &fs)
    {
        for (int phaseIdx = 0; phaseIdx < FluidState::numPhases; ++phaseIdx)
            updatePhase(fs, phaseIdx);
    };

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     */
    template <class FluidState>
    void updatePhase(const FluidState &fs, int phaseIdx)
    {};

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on temperature
     *
     * *Only* use this method if only the temperature of a phase
     * changed between two update*() calls. If more changed, call
     * updatePhase()!
     */
    template <class FluidState> void
    updatePhaseTemperature(const FluidState &fs, int phaseIdx)
    {
        asImp_().updatePhase(fs, phaseIdx);
    };

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on pressure
     *
     * *Only* use this method if only the pressure of a phase changed
     * between two update*() calls. If more changed, call
     * updatePhase()!
     */
    template <class FluidState>
    void updatePhasePressure(const FluidState &fs, int phaseIdx)
    {
        asImp_().updatePhase(fs, phaseIdx);
    };

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on composition
     *
     * *Only* use this method if neither the pressure nor the
     * temperature of the phase changed between two update*()
     * calls. If more changed, call updatePhase()!
     */
    template <class FluidState>
    void updatePhaseComposition(const FluidState &fs, int phaseIdx)
    {
        asImp_().updatePhase(fs, phaseIdx);
    };

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on the mole fraction of a single component
     *
     * *Only* use this method if just a single component's
     * concentration changed between two update*() calls. If more than
     * one concentration changed, call updatePhaseComposition() of
     * updatePhase()!
     */
    template <class FluidState>
    void updatePhaseSingleMoleFraction(const FluidState &fs,
                                       int phaseIdx,
                                       int compIdx)
    {
        asImp_().updatePhaseComposition(fs, phaseIdx);
    };

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
};

} // end namepace

#endif
