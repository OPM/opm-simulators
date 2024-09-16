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

#include <config.h>
#include <opm/models/nonlinear/newtonmethodparams.hpp>

#include <opm/models/utils/parametersystem.hpp>

#if HAVE_QUAD
#include <opm/material/common/quad.hpp>
#endif

namespace Opm {

template<class Scalar>
void NewtonMethodParams<Scalar>::registerParameters()
{
    Parameters::Register<Parameters::NewtonVerbose>
        ("Specify whether the Newton method should inform "
         "the user about its progress or not");
    Parameters::Register<Parameters::NewtonWriteConvergence>
        ("Write the convergence behaviour of the Newton "
         "method to a VTK file");
    Parameters::Register<Parameters::NewtonTargetIterations>
        ("The 'optimum' number of Newton iterations per time step");
    Parameters::Register<Parameters::NewtonMaxIterations>
        ("The maximum number of Newton iterations per time step");
    Parameters::Register<Parameters::NewtonTolerance<Scalar>>
        ("The maximum raw error tolerated by the Newton"
         "method for considering a solution to be converged");
    Parameters::Register<Parameters::NewtonMaxError<Scalar>>
        ("The maximum error tolerated by the Newton "
         "method to which does not cause an abort");
}

template<class Scalar>
void NewtonMethodParams<Scalar>::read()
{
    verbose_ = Parameters::Get<Parameters::NewtonVerbose>();
    writeConvergence_ = Parameters::Get<Parameters::NewtonWriteConvergence>();
    targetIterations_ = Parameters::Get<Parameters::NewtonTargetIterations>();
    maxIterations_ = Parameters::Get<Parameters::NewtonMaxIterations>();
    tolerance_ = Parameters::Get<Parameters::NewtonTolerance<Scalar>>();
    maxError_ = Parameters::Get<Parameters::NewtonMaxError<Scalar>>();
}

template struct NewtonMethodParams<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct NewtonMethodParams<float>;
#endif

#if HAVE_QUAD
template struct NewtonMethodParams<quad>;
#endif

} // namespace Opm
