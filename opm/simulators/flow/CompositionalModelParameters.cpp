/*
  Copyright 2026 NORCE Research AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/simulators/flow/CompositionalModelParameters.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

template <class Scalar>
void CompositionalModelParameters<Scalar>::registerParameters()
{
    BlackoilModelParameters<Scalar>::registerParameters();

    // Max. solution change tolerances
    Parameters::SetDefault<Parameters::ToleranceMaxDp<Scalar>>(1e4);
    Parameters::SetDefault<Parameters::ToleranceMaxDs<Scalar>>(1e-2);

    // Residual max-norm (CNV) and sum (MB) tolerances
    Parameters::SetDefault<Parameters::ToleranceCnv<Scalar>>(1e-4);
    Parameters::SetDefault<Parameters::ToleranceMb<Scalar>>(1e-20);
}

template struct CompositionalModelParameters<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct CompositionalModelParameters<float>;
#endif

} // namespace Opm
