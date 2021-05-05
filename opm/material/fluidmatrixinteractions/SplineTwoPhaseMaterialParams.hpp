// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::SplineTwoPhaseMaterialParams
 */
#ifndef OPM_SPLINE_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_SPLINE_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <opm/material/common/Spline.hpp>
#include <opm/material/common/EnsureFinalized.hpp>

#include <vector>
#include <cassert>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for a two-phase material law which
 *        uses a table and spline-based interpolation.
 */
template<class TraitsT>
class SplineTwoPhaseMaterialParams : public EnsureFinalized
{
    typedef typename TraitsT::Scalar Scalar;
public:
    using EnsureFinalized :: finalize;

public:
    typedef std::vector<Scalar> SamplePoints;
    typedef ::Opm::Spline<Scalar> Spline;
    typedef typename Spline::SplineType SplineType;

    typedef TraitsT Traits;

    SplineTwoPhaseMaterialParams()
    {
    }

    /*!
     * \brief Return the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const Spline& pcnwSpline() const
    { EnsureFinalized::check(); return pcwnSpline_; }

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
    { EnsureFinalized::check(); return krwSpline_; }

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
    { EnsureFinalized::check(); return krnSpline_; }

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
    Spline SwSpline_;
    Spline pcwnSpline_;
    Spline krwSpline_;
    Spline krnSpline_;
};
} // namespace Opm

#endif
