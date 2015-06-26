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

#ifndef OPM_BLACKOIL_STATE_TO_FLUID_STATE_HEADER_INCLUDED
#define OPM_BLACKOIL_STATE_TO_FLUID_STATE_HEADER_INCLUDED

#include <opm/core/simulator/BlackoilState.hpp>

namespace Opm
{

/*!
 * \brief This is an light weight "impedance adaption" class with a well defined API for
 *        saturation and PVT functions.
 *
 * It uses a stripped down version of opm-material's FluidState API and takes an
 * Opm::BlackoilState plus a cell index.
 *
 * Note that this class requires that is underlying BlackoilState must valid for at least
 * as long as an object of BlackoilStateToFluidState is used.
 */
class BlackoilStateToFluidState
{
public:
    typedef double Scalar;

    enum { numPhases = 3 };
    enum { numComponents = 3 };

    /*!
     * \brief Create a BlackoilState to Fluid state wrapper object.
     *
     * Note that this class requires that is underlying BlackoilState must valid for at least
     * as long as an object of BlackoilStateToFluidState is used.
     */
    BlackoilStateToFluidState(const BlackoilState& blackoilState)
        : blackoilState_(blackoilState)
    {
        if (blackoilState_.numPhases() != numPhases) {
            OPM_THROW(std::runtime_error,
                      "Only " << numPhases << " are supported, but the deck specifies " << blackoilState_.numPhases());
        }
    }

    /*!
     * \brief Sets the index of the currently used cell.
     *
     * After calling this, the values returned by the other methods are specific for this
     * cell.
     */
    void setCurrentCellIndex(unsigned cellIdx)
    { cellIdx_ = cellIdx; }

    /*!
     * \brief Returns the saturation of a phase for the current cell index.
     */
    Scalar saturation(int phaseIdx) const
    { return blackoilState_.saturation()[numPhases*cellIdx_ + phaseIdx]; }

    /*!
     * \brief Returns the temperature [K] of a phase for the current cell index.
     */
    Scalar temperature(int phaseIdx) const
    { return blackoilState_.temperature()[cellIdx_]; }

    // TODO (?) pressure, composition, etc

private:
    size_t numCells() const
    { return blackoilState_.pressure().size(); }

    const BlackoilState& blackoilState_;
    unsigned cellIdx_;
}

} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
