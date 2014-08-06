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
 * \copydoc Opm::SplineTwoPhaseMaterialParams
 */
#ifndef OPM_SPLINE_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_SPLINE_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <opm/core/utility/Spline.hpp>

#include <vector>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for a two-phase material law which
 *        uses a table and spline-based interpolation.
 */
template<class TraitsT>
class SplineTwoPhaseMaterialParams
{
    typedef typename TraitsT::Scalar Scalar;

public:
    typedef std::vector<Scalar> SamplePoints;
    typedef Opm::Spline<Scalar> Spline;
    typedef typename Spline::SplineType SplineType;

    typedef TraitsT Traits;

    SplineTwoPhaseMaterialParams()
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
    }

    /*!
     * \brief Return the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const Spline& pcnwSpline() const
    { assertFinalized_(); return pcwnSpline_; }

    /*!
     * \brief Set the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setPcnwSamples(const SamplePoints& SwSamplePoints,
                        const SamplePoints& pcnwSamplePoints,
                        SplineType splineType = Spline::Monotonic)
    {
        assert(SwSamplePoints.size() == pcnwSamplePoints.size());
        pcwnSpline_.setXYContainers(SwSamplePoints, pcnwSamplePoints, splineType);
    }

    /*!
     * \brief Return the sampling points for the relative permeability
     *        curve of the wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const Spline& krwSpline() const
    { assertFinalized_(); return krwSpline_; }

    /*!
     * \brief Set the sampling points for the relative permeability
     *        curve of the wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setKrwSamples(const SamplePoints& SwSamplePoints,
                       const SamplePoints& krwSamplePoints,
                       SplineType splineType = Spline::Monotonic)
    {
        assert(SwSamplePoints.size() == krwSamplePoints.size());
        krwSpline_.setXYContainers(SwSamplePoints, krwSamplePoints, splineType);
    }

    /*!
     * \brief Return the sampling points for the relative permeability
     *        curve of the non-wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const Spline& krnSpline() const
    { assertFinalized_(); return krnSpline_; }

    /*!
     * \brief Set the sampling points for the relative permeability
     *        curve of the non-wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setKrnSamples(const SamplePoints& SwSamplePoints,
                       const SamplePoints& krnSamplePoints,
                       SplineType splineType = Spline::Monotonic)
    {
        assert(SwSamplePoints.size() == krnSamplePoints.size());
        krnSpline_.setXYContainers(SwSamplePoints, krnSamplePoints, splineType);
    }

private:
#ifndef NDEBUG
    void assertFinalized_() const
    { assert(finalized_); }

    bool finalized_;
#else
    void assertFinalized_() const
    { }
#endif

    Spline SwSpline_;
    Spline pcwnSpline_;
    Spline krwSpline_;
    Spline krnSpline_;
};
} // namespace Opm

#endif
