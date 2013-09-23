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
 *
 * \brief Test for the isothermal VCVF discretization based on non-linear
 *        complementarity problems.
 */
#include "config.h"

#include <ewoms/common/start.hh>
#include <ewoms/models/ncp/ncpmodel.hh>
#include "problems/co2injectionproblem.hh"

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(Co2InjectionNcpProblem, INHERITS_FROM(VcfvNcp, Co2InjectionBaseProblem));
}
}

int main(int argc, char** argv)
{
    typedef TTAG(Co2InjectionNcpProblem) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
