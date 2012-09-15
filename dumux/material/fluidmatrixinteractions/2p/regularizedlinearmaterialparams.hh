// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          the linear constitutive relations.
 */
#ifndef REGULARIZED_LINEAR_PARAMS_HH
#define REGULARIZED_LINEAR_PARAMS_HH

#include "linearmaterialparams.hh"

namespace Dumux {
/*!
 * \ingroup fluidmatrixinteractionsparams
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          the linear constitutive relations.
 *
 */
template<class ScalarT>
class RegularizedLinearMaterialParams : public LinearMaterialParams<ScalarT>
{
public:
    typedef ScalarT Scalar;

    RegularizedLinearMaterialParams()
    {}

    /*!
     * \brief Return the threshold saturation respective phase below
     *        which the relative permeability gets regularized.
     *
     * This is just 5%. If you need a different value, write your own
     * parameter class.
     */
    Scalar krLowS() const
    { return 0.05; }

    /*!
     * \brief Return the threshold saturation of the respective phase
     *        above which the relative permeability gets regularized.
     *
     * This is just 95%. If you need a different value, write your own
     * parameter class.
     */
    Scalar krHighS() const
    { return 0.95; }

};
} // namespace Dumux

#endif
