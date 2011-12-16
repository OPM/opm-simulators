// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * Reference implementation of parameters for the M-phase brookscorey
 * material material.
 */
#ifndef DUMUX_MP_BROOKS_COREY_MATERIAL_PARAMS_HH
#define DUMUX_MP_BROOKS_COREY_MATERIAL_PARAMS_HH

namespace Dumux
{
/*!
 * \brief Reference implementation of a parameter class for the
 *        M-phase Brooks-Corey material law.
 */
template<int numPhasesV, class ScalarT>
class MpBrooksCoreyMaterialParams
{
public:
    typedef ScalarT Scalar;
    enum { numPhases = numPhasesV };

    MpBrooksCoreyMaterialParams()
    {
        for (int i = 0; i < numPhases; ++i) {
            setEntryPressure(i, 0.0);
            setLambda(i, 0.0);
        }
    }

    /*!
     * \brief Return the entry pressure for a phase.
     */
    Scalar entryPressure(int phaseIdx) const
    { return entryPressure_[phaseIdx]; }

    /*!
     * \brief Set the entry pressure for a phase.
     */
    void setEntryPressure(int phaseIdx, Scalar val) const
    { entryPressure_[phaseIdx] = val; }

    /*!
     * \brief Return the alpha shape parameter of the Brooks-Corey law
     *        for a phase.
     */
    Scalar alpha(int phaseIdx) const
    { return alpha_[phaseIdx]; }

    /*!
     * \brief Set the alpha shape parameter of the Brooks-Corey law
     *        for a phase.
     */
    void setLambda(int phaseIdx, Scalar val) const
    { alpha_[phaseIdx] = val; }

private:
    Scalar entryPressure_[numPhases];
    Scalar alpha_[numPhases];
};
} // namespace Dumux

#endif
