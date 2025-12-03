// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef ELASTICITY_INDICES_HPP
#define ELASTICITY_INDICES_HPP


namespace Opm {

template <unsigned PVOffset>
struct ElasticityIndices
{
    // ///
    // Primary variable indices
    // ///
    // Starting index of displacement vector
    static constexpr int disp0Idx = PVOffset + 0;

    // Starting index of rotation vector
    static constexpr int rot0Idx = disp0Idx + 3;

    // Index of solid pressure
    static constexpr int solidPres0Idx = rot0Idx + 3;

    // ///
    // Equation indices
    // ///
    // Starting index of all equations, in particular the stress equations
    static constexpr int conti0EqIdx = PVOffset + 0;

    // Starting index of rotation equations
    static constexpr int contiRotEqIdx = conti0EqIdx + 3;

    // Starting index of solid pressure equation
    static constexpr int contiSolidPresEqIdx = contiRotEqIdx + 3;

    // Total number of equations/primary variables
    static constexpr int numEq = 7;
};  // struct ElasticityIndices

}  // namespace Opm

#endif
