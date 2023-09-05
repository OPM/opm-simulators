/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <string>
#include <vector>

namespace Opm {
class ConvergenceReport;
}

namespace Opm {
namespace detail {

template<class Scalar>
ConvergenceReport createConvergenceReport(const int iteration,
                                          const double reportTime,
                                          const double tol_mb,
                                          const double tol_cnv,
                                          const double maxResidualAllowed,
                                          const bool terminal_output,
                                          const std::vector<std::string>& names,
                                          const std::vector<Scalar>& CNV,
                                          const std::vector<Scalar>& mass_balance_residual,
                                          const std::string& terminal_header);

template<class Scalar>
std::tuple<double,double> convergenceReduction(Parallel::Communication comm,
                                               const double pvSumLocal,
                                               const double numAquiferPvSumLocal,
                                               std::vector<Scalar>& R_sum,
                                               std::vector<Scalar>& maxCoeff,
                                               std::vector<Scalar>& B_avg);

} // namespace detail
} // namespace Opm
