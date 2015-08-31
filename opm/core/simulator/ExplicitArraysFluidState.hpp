/*
  Copyright 2015 Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_EXPLICIT_ARRAYS_FLUID_STATE_HEADER_INCLUDED
#define OPM_EXPLICIT_ARRAYS_FLUID_STATE_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>

#include <array>

namespace Opm
{

/*!
 * \brief This is a fluid state which translates global arrays and translates them to a
 *        subset of the fluid state API.
 *
 * This class is similar to Opm::BlackoilStateToFluidState.
 */
class ExplicitArraysFluidState
{
public:
    typedef double Scalar;
    enum { numPhases = BlackoilPhases::MaxNumPhases };

    explicit ExplicitArraysFluidState(const PhaseUsage& phaseUsage)
        : phaseUsage_(phaseUsage)
    {}

    /*!
     * \brief Sets the currently used array index.
     *
     * After calling this, the values returned by the other methods are specific for this
     * index.
     */
    void setIndex(unsigned arrayIdx)
    {
        int np = phaseUsage_.num_phases;
        for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
            if (!phaseUsage_.phase_used[phaseIdx]) {
                sats_[phaseIdx] = 0.0;
            }
            else {
                sats_[phaseIdx] = saturations_[np*arrayIdx + phaseUsage_.phase_pos[phaseIdx]];
            }
        }
    }

    /*!
     * \brief Set the array containing the phase saturations.
     *
     * This array is supposed to be of size numPhases*size and is not allowed to be
     * deleted while the ExplicitArraysFluidState object is alive. This class assumes
     * that the saturations of all phase saturations for a point are consequtive, i.e.,
     * in the array the saturations cycle fastest.
     */
    void setSaturationArray(const double* saturations)
    { saturations_ = saturations; }

    /*!
     * \brief Returns the saturation of a phase for the current cell index.
     */
    Scalar saturation(int phaseIdx) const
    { return sats_[phaseIdx]; }

    // TODO (?) temperature, pressure, composition, etc

private:
    const PhaseUsage phaseUsage_;
    const double* saturations_;
    std::array<Scalar, BlackoilPhases::MaxNumPhases> sats_;
};

} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
