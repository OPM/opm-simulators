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

#include <opm/core/simulator/BlackoilState.hpp>

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

    explicit ExplicitArraysFluidState(const unsigned int num_phases)
	: numPhases_(num_phases)
    {}

    /*!
     * \brief Sets the currently used array index.
     *
     * After calling this, the values returned by the other methods are specific for this
     * index.
     */
    void setIndex(unsigned arrayIdx)
    { arrayIdx_ = arrayIdx; }

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
     * \brief Set the array containing the phase temperatures.
     *
     * This array is supposed to be of size 'size' and is not allowed to be
     * deleted while the ExplicitArraysFluidState object is alive.
     */
    void setTemperatureArray(const double* temperature)
    { temperature_ = temperature; }

    /*!
     * \brief Returns the saturation of a phase for the current cell index.
     */
    Scalar saturation(int phaseIdx) const
    { return saturations_[numPhases_*arrayIdx_ + phaseIdx]; }

    /*!
     * \brief Returns the temperature [K] of a phase for the current cell index.
     */
    Scalar temperature(int /* phaseIdx */) const
    { return temperature_[arrayIdx_]; }

    // TODO (?) pressure, composition, etc

private:
    const double* saturations_;
    const double* temperature_;

    unsigned arrayIdx_;
    unsigned numPhases_;
};

} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
