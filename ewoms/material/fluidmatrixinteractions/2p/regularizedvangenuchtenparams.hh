// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Philipp Nuske                                     *
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
 * \copydoc Ewoms::RegularizedVanGenuchtenParams
 */
#ifndef EWOMS_REGULARIZED_VAN_GENUCHTEN_PARAMS_HH
#define EWOMS_REGULARIZED_VAN_GENUCHTEN_PARAMS_HH

#include "vangenuchtenparams.hh"

namespace Ewoms {
/*!
 * \ingroup fluidmatrixinteractionsparams
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          VanGenuchten "material law".
 *
 */
template<class ScalarT>
class RegularizedVanGenuchtenParams : public VanGenuchtenParams<ScalarT>
{
public:
    typedef ScalarT Scalar;
    typedef VanGenuchtenParams<Scalar> Parent;

    RegularizedVanGenuchtenParams()
        : pCLowSw_(0.01)
        , pCHighSw_(0.99)
    {}

    RegularizedVanGenuchtenParams(Scalar vgAlpha,
                                   Scalar vgN)
        : Parent(vgAlpha, vgN)
        , pCLowSw_(0.01)
        , pCHighSw_(0.99)
    {}

    /*!
     * \brief Return the threshold saturation below which the
     *        capillary pressure is regularized.
     */
    Scalar pCLowSw() const
    { return pCLowSw_; }

    /*!
     * \brief Set the threshold saturation below which the capillary
     *        pressure is regularized.
     */
    void setPCLowSw(Scalar value)
    { pCLowSw_ = value; }

    /*!
     * \brief Return the threshold saturation below which the
     *        capillary pressure is regularized.
     */
    Scalar pCHighSw() const
    { return pCHighSw_; }

    /*!
     * \brief Set the threshold saturation below which the capillary
     *        pressure is regularized.
     */
    void setPCHighSw(Scalar value)
    { pCHighSw_ = value; }

    /*!
     * \brief Threshold saturation below which the relative
     *        permeability of the non-wetting phase gets regulatized.
     *
     * This is just 10%. If you need a different value, overload this
     * class.
     */
    Scalar krnLowSw() const
    { return 0.10; }

    /*!
     * \brief Threshold saturation above which the relative
     *        permeability of the wetting phase gets regulatized.
     *
     * This is just 90%. If you need a different value, overload this
     * class.
     */
    Scalar krwHighSw() const
    { return 0.90; }

private:
    Scalar pCLowSw_;
    Scalar pCHighSw_;
};
} // namespace Ewoms

#endif
