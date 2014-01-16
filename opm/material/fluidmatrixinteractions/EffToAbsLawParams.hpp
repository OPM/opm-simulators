/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 * \copydoc Opm::EffToAbsLawParams
 */
#ifndef OPM_EFF_TO_ABS_LAW_PARAMS_HPP
#define OPM_EFF_TO_ABS_LAW_PARAMS_HPP

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws from effective to absolute
 *        saturations.
 */
template <class EffLawParamsT, int numPhases>
class EffToAbsLawParams : public EffLawParamsT
{
    typedef EffLawParamsT EffLawParams;
    typedef typename EffLawParams::Traits::Scalar Scalar;

public:
    typedef typename EffLawParams::Traits Traits;


    EffToAbsLawParams()
        : EffLawParams()
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            residualSaturation_[phaseIdx] = 0.0;

#ifndef NDEBUG
        finalized_ = false;
#endif
    }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        sumResidualSaturations_ = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            sumResidualSaturations_ += residualSaturation_[phaseIdx];

        EffLawParams::finalize();

#ifndef NDEBUG
        finalized_ = true;
#endif
    }

    /*!
     * \brief Return the residual saturation of a phase.
     */
    Scalar residualSaturation(int phaseIdx) const
    { assertFinalized_(); return residualSaturation_[phaseIdx]; }

    /*!
     * \brief Return the sum of the residual saturations.
     */
    Scalar sumResidualSaturations() const
    { assertFinalized_(); return sumResidualSaturations_; }

    /*!
     * \brief Set the residual saturation of a phase.
     */
    void setResidualSaturation(int phaseIdx, Scalar value)
    { residualSaturation_[phaseIdx] = value; }

private:
#ifndef NDEBUG
    void assertFinalized_() const
    { assert(finalized_); }

    bool finalized_;
#else
    void assertFinalized_() const
    { }
#endif

    Scalar residualSaturation_[numPhases];
    Scalar sumResidualSaturations_;
};

} // namespace Opm

#endif
