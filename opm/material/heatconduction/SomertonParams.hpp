/*
  Copyright (C) 2012-2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Opm::SomertonParams
 */
#ifndef OPM_SOMERTON_PARAMS_HPP
#define OPM_SOMERTON_PARAMS_HPP

#include <cassert>

namespace Opm {

/*!
 * \brief The default implementation of a parameter object for the
 *        Somerton heatconduction law.
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

} // namespace Opm

#endif
