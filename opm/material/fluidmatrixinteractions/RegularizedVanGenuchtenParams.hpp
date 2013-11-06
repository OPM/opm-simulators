// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 * \copydoc Opm::RegularizedVanGenuchtenParams
 */
#ifndef OPM_REGULARIZED_VAN_GENUCHTEN_PARAMS_HH
#define OPM_REGULARIZED_VAN_GENUCHTEN_PARAMS_HH

#include "VanGenuchtenParams.hpp"

#include <opm/core/utility/Spline.hpp>

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
    typedef Opm::VanGenuchten<TraitsT> VanGenuchten;

public:
    typedef TraitsT Traits;

    RegularizedVanGenuchtenParams()
        : pcwnLowSw_(0.01)
        , pcwnHighSw_(0.99)
    {}

    RegularizedVanGenuchtenParams(Scalar vgAlpha, Scalar vgN)
        : Parent(vgAlpha, vgN)
        , pcwnLowSw_(0.01)
        , pcwnHighSw_(0.99)
    {}

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        Parent::finalize();

        pcwnLow_ = VanGenuchten::twoPhaseSatPcwn(*this, pcwnLowSw_);
        pcwnSlopeLow_ = VanGenuchten::twoPhaseSatDpcwn_dSw(*this, pcwnLowSw_);
        pcwnHigh_ = VanGenuchten::twoPhaseSatPcwn(*this, pcwnHighSw_);
        pcwnSlopeHigh_ = 2*(0.0 - pcwnHigh_)/(1.0 - pcwnHighSw_);

        Scalar mThreshold = VanGenuchten::twoPhaseSatDpcwn_dSw(*this, pcwnHighSw_);

        pcwnHighSpline_.set(pcwnHighSw_, 1.0, // x0, x1
                            pcwnHigh_, 0, // y0, y1
                            mThreshold, pcwnSlopeHigh_); // m0, m1

    }

    /*!
     * \brief Return the threshold saturation below which the
     *        capillary pressure is regularized.
     */
    Scalar pcwnLowSw() const
    { return pcwnLowSw_; }

    /*!
     * \brief Return the capillary pressure at the low threshold
     *        saturation of the wetting phase.
     */
    Scalar pcwnLow() const
    { return pcwnLow_; }

    /*!
     * \brief Return the slope capillary pressure curve if Sw is
     *        smaller or equal to the low threshold saturation.
     *
     * For this case, we extrapolate the curve using a straight line.
     */
    Scalar pcwnSlopeLow() const
    { return pcwnSlopeLow_; }

    /*!
     * \brief Set the threshold saturation below which the capillary
     *        pressure is regularized.
     */
    void setPCLowSw(Scalar value)
    { pcwnLowSw_ = value; }

    /*!
     * \brief Return the threshold saturation below which the
     *        capillary pressure is regularized.
     */
    Scalar pcwnHighSw() const
    { return pcwnHighSw_; }

    /*!
     * \brief Return the capillary pressure at the high threshold
     *        saturation of the wetting phase.
     */
    Scalar pcwnHigh() const
    { return pcwnHigh_; }

    /*!
     * \brief Return the spline curve which ought to be used between
     *        the upper threshold saturation and 1.
     */
    const Spline<Scalar> &pcwnHighSpline() const
    { return pcwnHighSpline_; }

    /*!
     * \brief Return the slope capillary pressure curve if Sw is
     *        larger or equal to 1.
     *
     * For this case, we extrapolate the curve using a straight line.
     */
    Scalar pcwnSlopeHigh() const
    { return pcwnSlopeHigh_; }

    /*!
     * \brief Set the threshold saturation below which the capillary
     *        pressure is regularized.
     */
    void setPCHighSw(Scalar value)
    { pcwnHighSw_ = value; }

private:
    Scalar pcwnLowSw_;
    Scalar pcwnHighSw_;

    Scalar pcwnLow_;
    Scalar pcwnHigh_;

    Scalar pcwnSlopeLow_;
    Scalar pcwnSlopeHigh_;

    Spline<Scalar> pcwnHighSpline_;
};
} // namespace Opm

#endif
