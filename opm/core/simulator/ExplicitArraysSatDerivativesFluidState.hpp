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

#ifndef OPM_EXPLICIT_ARRAYS_SAT_DERIVATIVES_FLUID_STATE_HEADER_INCLUDED
#define OPM_EXPLICIT_ARRAYS_SAT_DERIVATIVES_FLUID_STATE_HEADER_INCLUDED

#include <opm/material/localad/Evaluation.hpp>
#include <opm/material/localad/Math.hpp>

namespace Opm
{

// this class does not need to be implemented. its only purpose is to make compiler
// messages for local-AD framework more explicit (because the class name of the
// Evaluation will include a tag name so that it is clear which derivatives are handled).
class SaturationDerivativesTag
{};

/*!
 * \brief This is a fluid state which translates global arrays and translates them to a
 *        subset of the fluid state API.
 *
 * This class is similar to Opm::ExplicitArraysFluidState except for the fact that it
 * uses opm-material's local automatic differentiation framework to allow implicit
 * treatment of saturation derivatives.
 */
class ExplicitArraysSatDerivativesFluidState
{
public:
    enum { numPhases = BlackoilPhases::MaxNumPhases };
    enum { numComponents = 3 };

    typedef Opm::LocalAd::Evaluation<double, SaturationDerivativesTag, numPhases> Evaluation;
    typedef Evaluation Scalar;

    ExplicitArraysSatDerivativesFluidState(const PhaseUsage& phaseUsage)
        : phaseUsage_(phaseUsage)
    {
        globalSaturationArray_ = 0;

        // initialize the evaluation objects for the saturations
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            saturation_[phaseIdx] = Evaluation::createVariable(0.0, phaseIdx);
        }
    }

    /*!
     * \brief Sets the currently used array index.
     *
     * After calling this, the values returned by the other methods are specific for this
     * index.
     */
    void setIndex(unsigned arrayIdx)
    {
        int np = phaseUsage_.num_phases;

        // copy the saturations values from the global value. the derivatives do not need
        // to be modified for these...
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!phaseUsage_.phase_used[phaseIdx]) {
                saturation_[phaseIdx].value = 0.0;
            }
            else {
                saturation_[phaseIdx].value = globalSaturationArray_[np*arrayIdx + phaseUsage_.phase_pos[phaseIdx]];
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
    { globalSaturationArray_ = saturations; }

    /*!
     * \brief Returns the saturation of a phase for the current cell index.
     */
    const Evaluation& saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    // TODO (?) temperature, pressure, composition, etc

private:
    const PhaseUsage phaseUsage_;
    const double* globalSaturationArray_;

    std::array<Evaluation, numPhases> saturation_;
};

} // namespace Opm

#endif
