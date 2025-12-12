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
#include "config.h"

#include <opm/models/tpsa/tpsanewtonmethodparams.hpp>
#include <opm/models/utils/parametersystem.hpp>

#if HAVE_QUAD
#include <opm/material/common/quad.hpp>
#endif


namespace Opm {

/*!
* \brief Register runtime parameters for TPSA Newton method
*/
template<class Scalar>
void TpsaNewtonMethodParams<Scalar>::registerParameters()
{
    Parameters::Register<Parameters::TpsaNewtonVerbose>
        ("Specify whether the TPSA Newton method should inform the user about its progress or not");
    Parameters::Register<Parameters::TpsaNewtonWriteConvergence>
        ("Write the convergence behaviour of the TPSA Newton method to a VTK file");
    Parameters::Register<Parameters::TpsaNewtonTargetIterations>
        ("The 'optimum' number of TPSA Newton iterations");
    Parameters::Register<Parameters::TpsaNewtonMaxIterations>
        ("The maximum number of TPSA Newton iterations");
    Parameters::Register<Parameters::TpsaNewtonMinIterations>
        ("The minimum number of TPSA Newton iterations");
    Parameters::Register<Parameters::TpsaNewtonTolerance<Scalar>>
        ("The maximum raw error tolerated by the TPSA Newton method for considering a solution to be converged");
    Parameters::Register<Parameters::TpsaNewtonMaxError<Scalar>>
        ("The maximum error tolerated by the TPSA Newton method to which does not cause an abort");
}

/*!
* \brief Read and internalize runtime parameters for TPSA Newton method
*/
template<class Scalar>
void TpsaNewtonMethodParams<Scalar>::read()
{
    verbose_ = Parameters::Get<Parameters::TpsaNewtonVerbose>();
    writeConvergence_ = Parameters::Get<Parameters::TpsaNewtonWriteConvergence>();
    targetIterations_ = Parameters::Get<Parameters::TpsaNewtonTargetIterations>();
    minIterations_ = Parameters::Get<Parameters::TpsaNewtonMinIterations>();
    maxIterations_ = Parameters::Get<Parameters::TpsaNewtonMaxIterations>();
    tolerance_ = Parameters::Get<Parameters::TpsaNewtonTolerance<Scalar>>();
    maxError_ = Parameters::Get<Parameters::TpsaNewtonMaxError<Scalar>>();
}

template struct TpsaNewtonMethodParams<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct TpsaNewtonMethodParams<float>;
#endif

#if HAVE_QUAD
template struct TpsaNewtonMethodParams<quad>;
#endif

} // namespace Opm
