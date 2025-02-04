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
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_MECH_CONTAINER_HPP
#define OPM_MECH_CONTAINER_HPP

#include <vector>

namespace Opm {

template<class Scalar>
class MechContainer
{
    using ScalarBuffer = std::vector<Scalar>;

public:
    ScalarBuffer potentialForce_;
    ScalarBuffer potentialPressForce_;
    ScalarBuffer potentialTempForce_;

    ScalarBuffer dispX_;
    ScalarBuffer dispY_;
    ScalarBuffer dispZ_;
    ScalarBuffer stressXX_;
    ScalarBuffer stressYY_;
    ScalarBuffer stressZZ_;
    ScalarBuffer stressXY_;
    ScalarBuffer stressXZ_;
    ScalarBuffer stressYZ_;
    ScalarBuffer delstressXX_;
    ScalarBuffer delstressYY_;
    ScalarBuffer delstressZZ_;
    ScalarBuffer delstressXY_;
    ScalarBuffer delstressXZ_;
    ScalarBuffer delstressYZ_;
    ScalarBuffer strainXX_;
    ScalarBuffer strainYY_;
    ScalarBuffer strainZZ_;
    ScalarBuffer strainXY_;
    ScalarBuffer strainXZ_;
    ScalarBuffer strainYZ_;
    ScalarBuffer linstressXX_;
    ScalarBuffer linstressYY_;
    ScalarBuffer linstressZZ_;
    ScalarBuffer linstressXY_;
    ScalarBuffer linstressXZ_;
    ScalarBuffer linstressYZ_;
    ScalarBuffer fracstressXX_;
    ScalarBuffer fracstressYY_;
    ScalarBuffer fracstressZZ_;
    ScalarBuffer fracstressXY_;
    ScalarBuffer fracstressXZ_;
    ScalarBuffer fracstressYZ_;
};

} // namespace Opm

#endif // OPM_MECH_CONTAINER_HPP
