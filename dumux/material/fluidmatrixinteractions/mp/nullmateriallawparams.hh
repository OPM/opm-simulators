// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \file linearmaterialparams.hh
 *
 * Reference implementation of parameters for the M-phase linear
 * material material.
 */
#ifndef DUMUX_NULL_MATERIAL_LAW_PARAMS_HH
#define DUMUX_NULL_MATERIAL_LAW_PARAMS_HH

namespace Dumux {
/*!
 * \brief Reference implementation of params for the linear M-phase
 *        material material.
 */
template<int numPhasesV, class ScalarT>
class NullMaterialLawParams
{
public:
    typedef ScalarT Scalar;
    enum { numPhases = numPhasesV };

    NullMaterialLawParams()
    { }
};
} // namespace Dumux

#endif
