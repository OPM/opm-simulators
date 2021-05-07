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
 * \copydoc Opm::RegularizedVanGenuchtenParams
 */
#ifndef OPM_REGULARIZED_VAN_GENUCHTEN_PARAMS_HPP
#define OPM_REGULARIZED_VAN_GENUCHTEN_PARAMS_HPP

#include "VanGenuchten.hpp"
#include "VanGenuchtenParams.hpp"

#include <opm/material/common/Spline.hpp>

#include <cassert>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          VanGenuchten "material law".
 *
 */
template<class TraitsT>
class RegularizedVanGenuchtenParams : public VanGenuchtenParams<TraitsT>
{
    typedef typename TraitsT::Scalar Scalar;
    typedef VanGenuchtenParams<TraitsT> Parent;
    typedef ::Opm::VanGenuchten<TraitsT> VanGenuchten;

public:
    using Parent :: finalize;

    typedef TraitsT Traits;

    RegularizedVanGenuchtenParams()
        : pcnwLowSw_(0.01)
        , pcnwHighSw_(0.99)
    {}

    RegularizedVanGenuchtenParams(Scalar vgAlpha, Scalar vgN)
        : Parent(vgAlpha, vgN)
        , pcnwLowSw_(0.01)
        , pcnwHighSw_(0.99)
    {
        finalize();
    }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        Parent::finalize();

        pcnwLow_ = VanGenuchten::twoPhaseSatPcnw(*this, pcnwLowSw_);
        pcnwSlopeLow_ = dPcnw_dSw_(pcnwLowSw_);
        pcnwHigh_ = VanGenuchten::twoPhaseSatPcnw(*this, pcnwHighSw_);
        pcnwSlopeHigh_ = 2*(0.0 - pcnwHigh_)/(1.0 - pcnwHighSw_);

        Scalar mThreshold = dPcnw_dSw_(pcnwHighSw_);

        pcnwHighSpline_.set(pcnwHighSw_, 1.0, // x0, x1
                            pcnwHigh_, 0, // y0, y1
                            mThreshold, pcnwSlopeHigh_); // m0, m1
    }

    /*!
     * \brief Return the threshold saturation below which the
     *        capillary pressure is regularized.
     */
    Scalar pcnwLowSw() const
    { EnsureFinalized::check(); return pcnwLowSw_; }

    /*!
     * \brief Return the capillary pressure at the low threshold
     *        saturation of the wetting phase.
     */
    Scalar pcnwLow() const
    { EnsureFinalized::check(); return pcnwLow_; }

    /*!
     * \brief Return the slope capillary pressure curve if Sw is
     *        smaller or equal to the low threshold saturation.
     *
     * For this case, we extrapolate the curve using a straight line.
     */
    Scalar pcnwSlopeLow() const
    { EnsureFinalized::check(); return pcnwSlopeLow_; }

    /*!
     * \brief Set the threshold saturation below which the capillary
     *        pressure is regularized.
     */
    void setPCLowSw(Scalar value)
    { pcnwLowSw_ = value; }

    /*!
     * \brief Return the threshold saturation below which the
     *        capillary pressure is regularized.
     */
    Scalar pcnwHighSw() const
    { EnsureFinalized::check(); return pcnwHighSw_; }

    /*!
     * \brief Return the capillary pressure at the high threshold
     *        saturation of the wetting phase.
     */
    Scalar pcnwHigh() const
    { EnsureFinalized::check(); return pcnwHigh_; }

    /*!
     * \brief Return the spline curve which ought to be used between
     *        the upper threshold saturation and 1.
     */
    const Spline<Scalar>& pcnwHighSpline() const
    { EnsureFinalized::check(); return pcnwHighSpline_; }

    /*!
     * \brief Return the slope capillary pressure curve if Sw is
     *        larger or equal to 1.
     *
     * For this case, we extrapolate the curve using a straight line.
     */
    Scalar pcnwSlopeHigh() const
    { EnsureFinalized::check(); return pcnwSlopeHigh_; }

    /*!
     * \brief Set the threshold saturation below which the capillary
     *        pressure is regularized.
     */
    void setPCHighSw(Scalar value)
    { pcnwHighSw_ = value; }

private:
    Scalar dPcnw_dSw_(Scalar Sw) const
    {
        // use finite differences to calculate the derivative w.r.t. Sw of the
        // unregularized curve's capillary pressure.
        const Scalar eps = 1e-7;
        Scalar pc1 = VanGenuchten::twoPhaseSatPcnw(*this, Sw - eps);
        Scalar pc2 = VanGenuchten::twoPhaseSatPcnw(*this, Sw + eps);
        return (pc2 - pc1)/(2*eps);
    }

    Scalar pcnwLowSw_;
    Scalar pcnwHighSw_;

    Scalar pcnwLow_;
    Scalar pcnwHigh_;

    Scalar pcnwSlopeLow_;
    Scalar pcnwSlopeHigh_;

    Spline<Scalar> pcnwHighSpline_;
};
} // namespace Opm

#endif
