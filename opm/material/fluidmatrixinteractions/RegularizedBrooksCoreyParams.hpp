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
 * \copydoc Opm::RegularizedBrooksCoreyParams
 */
#ifndef OPM_REGULARIZED_BROOKS_COREY_PARAMS_HPP
#define OPM_REGULARIZED_BROOKS_COREY_PARAMS_HPP

#include "BrooksCoreyParams.hpp"

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          the Brooks-Corey capillary pressure model.
 */
template <class TraitsT>
class RegularizedBrooksCoreyParams : public Opm::BrooksCoreyParams<TraitsT>
{
    typedef Opm::BrooksCoreyParams<TraitsT> BrooksCoreyParams;
    typedef typename TraitsT::Scalar Scalar;

public:
    typedef TraitsT Traits;

    RegularizedBrooksCoreyParams()
        : BrooksCoreyParams()
        , SwThres_(1e-2)
    { }

    RegularizedBrooksCoreyParams(Scalar entryPressure, Scalar lambda)
        : BrooksCoreyParams(entryPressure, lambda)
        , SwThres_(1e-2)
    { }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        BrooksCoreyParams::finalize();
    }

    /*!
     * \brief Return the threshold saturation below which the capillary pressure
     *        is regularized.
     */
    Scalar thresholdSw() const
    { return SwThres_; }

    /*!
     * \brief Set the threshold saturation below which the capillary pressure
     *        is regularized.
     */
    void setThresholdSw(Scalar value)
    { SwThres_ = value; }

private:
    Scalar SwThres_;
};
} // namespace Opm

#endif
