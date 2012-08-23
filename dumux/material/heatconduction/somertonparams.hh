// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \file somertonparams.hh
 *
 * Reference implementation of parameters for the Somerton law of
 * heat conductivity in a porous medium.
 *
 * See:
 *
 * W.H. Somerton, A.H. El-Shaarani and S.M. Mobarak: High
 * Temperature Behavior of Rocks Associated with Geothermal Type
 * Reservoirs, paper SPE-4897 presentet at SPE California Regional
 * Meeting 1974, 1974
 *
 * or
 *
 * H. Class: Theorie und numerische Modellierung nichtisothermer
 * Mehrphasenprozesse in NAPL kontaminierten poroesen Medien, PhD
 * thesis, Technical University of Braunschweig, 2000
 */
#ifndef DUMUX_SOMERTON_PARAMS_HH
#define DUMUX_SOMERTON_PARAMS_HH

#include <cassert>

namespace Dumux {
/*!
 * \brief
 */
template <int numPhases, class ScalarT>
class SomertonParams
{
    // do not copy!
    SomertonParams(const SomertonParams &)
    {}

public:
    typedef ScalarT Scalar;

    SomertonParams()
    { }

    /*!
     * \brief Return the "fully saturated" heat conductivity of the
     *        porous medium [W/m^2 / (K/m)].
     *
     * In this context "fully saturated" means that the whole pore
     * space of the porous medium is filled by a given fluid phase.
     */
    Scalar fullySaturatedLambda(int phaseIdx) const
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return fullySaturatedLambda_[phaseIdx];
    }

    /*!
     * \brief Set the "fully saturated" heat conductivity of the
     *        porous medium [W/m^2 / (K/m)].
     *
     * In this context "fully saturated" means that the whole pore
     * space of the porous medium is filled by a given fluid phase.
     */
    void setFullySaturatedLambda(int phaseIdx, Scalar value)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(value > 0);

        fullySaturatedLambda_[phaseIdx] = value;
    }

    /*!
     * \brief Return the heat conductivity of the porous medium at
     *        vacuum [W/m^2 / (K/m)].
     */
    Scalar vacuumLambda() const
    {
        return vacuumLambda_;
    }

    /*!
     * \brief Set the "fully saturated" heat conductivity of the
     *        porous medium [W/m^2 / (K/m)].
     *
     * In this context "fully saturated" means that the whole pore
     * space of the porous medium is filled by a given fluid phase.
     */
    void setVacuumLambda(Scalar value)
    {
        assert(value > 0);

        vacuumLambda_ = value;
    }

private:
    Scalar fullySaturatedLambda_[numPhases];
    Scalar vacuumLambda_;
};

} // namespace Dumux

#endif
