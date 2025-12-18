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
#ifndef TPSA_NEWTON_METHOD_PARAMS_HPP
#define TPSA_NEWTON_METHOD_PARAMS_HPP


namespace Opm::Parameters {

// The maximum error which may occur in a simulation before the Newton method for the time step is aborted
template<class Scalar>
struct TpsaNewtonMaxError { static constexpr Scalar value = 1e100; };

// Number of maximum iterations for the Newton method
struct TpsaNewtonMaxIterations { static constexpr int value = 20; };

// Number of minimum iterations for the Newton method
struct TpsaNewtonMinIterations { static constexpr int value = 1; };

// Target number of iterations
struct TpsaNewtonTargetIterations { static constexpr int value = 10; };

// Convergence tolerance
template<class Scalar>
struct TpsaNewtonTolerance { static constexpr Scalar value = 1e-3; };

// Specifies verbosity of print messages
struct TpsaNewtonVerbosity { static constexpr int value = 1; };

} // namespace Opm::Parameters

namespace Opm {

/*!
* \brief Struct holding the parameters for TpsaNewtonMethod.
*/
template<class Scalar>
struct TpsaNewtonMethodParams
{
    static void registerParameters();
    void read();

    int verbosity_;
    bool writeConvergence_;
    int targetIterations_;
    int minIterations_;
    int maxIterations_;
    Scalar tolerance_;
    Scalar maxError_;
};  // struct TpsaNewtonMethodParams

}  // namespace Opm

#endif
