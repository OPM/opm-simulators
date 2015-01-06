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
 * \copydoc Opm::PiecewiseLinearTwoPhaseMaterialParams
 */
#ifndef OPM_PIECEWISE_LINEAR_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_PIECEWISE_LINEAR_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <vector>

#include <cassert>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for a two-phase material law which
 *        uses a table and piecewise constant interpolation.
 */
template<class TraitsT>
class PiecewiseLinearTwoPhaseMaterialParams
{
    typedef typename TraitsT::Scalar Scalar;

public:
    typedef std::vector<Scalar> ValueVector;

    typedef TraitsT Traits;

    PiecewiseLinearTwoPhaseMaterialParams()
    {
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
#ifndef NDEBUG
        finalized_ = true;
#endif

        // revert the order of the sampling points if they were given
        // in reverse direction.
        if (SwSamples_.front() > SwSamples_.back()) {
            for (unsigned origSampleIdx = 0;
                 origSampleIdx < SwSamples_.size() / 2;
                 ++ origSampleIdx)
            {
                unsigned newSampleIdx = SwSamples_.size() - origSampleIdx - 1;

                std::swap(SwSamples_[origSampleIdx], SwSamples_[newSampleIdx]);
                std::swap(pcwnSamples_[origSampleIdx], pcwnSamples_[newSampleIdx]);
                std::swap(krwSamples_[origSampleIdx], krwSamples_[newSampleIdx]);
                std::swap(krnSamples_[origSampleIdx], krnSamples_[newSampleIdx]);
            }
        }
    }

    /*!
     * \brief Return the wetting-phase saturation values of all sampling points.
     */
    const ValueVector& SwSamples() const
    { assertFinalized_(); return SwSamples_; }

    /*!
     * \brief Return the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const ValueVector& pcnwSamples() const
    { assertFinalized_(); return pcwnSamples_; }

    /*!
     * \brief Set the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setPcnwSamples(const ValueVector& SwValues, const ValueVector& values)
    { SwSamples_ = SwValues; pcwnSamples_ = values; }

    /*!
     * \brief Return the sampling points for the relative permeability
     *        curve of the wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const ValueVector& krwSamples() const
    { assertFinalized_(); return krwSamples_; }

    /*!
     * \brief Set the sampling points for the relative permeability
     *        curve of the wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setKrwSamples(const ValueVector& SwValues, const ValueVector& values)
    { SwSamples_ = SwValues; krwSamples_ = values; }

    /*!
     * \brief Return the sampling points for the relative permeability
     *        curve of the non-wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const ValueVector& krnSamples() const
    { assertFinalized_(); return krnSamples_; }

    /*!
     * \brief Set the sampling points for the relative permeability
     *        curve of the non-wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setKrnSamples(const ValueVector& SwValues, const ValueVector& values)
    { SwSamples_ = SwValues; krnSamples_ = values; }

private:
#ifndef NDEBUG
    void assertFinalized_() const
    { assert(finalized_); }

    bool finalized_;
#else
    void assertFinalized_() const
    { }
#endif

    ValueVector SwSamples_;
    ValueVector pcwnSamples_;
    ValueVector krwSamples_;
    ValueVector krnSamples_;
};
} // namespace Opm

#endif
