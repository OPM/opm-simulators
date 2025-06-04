/*
  Copyright 2025 Equinor ASA

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
#ifndef OPM_PRINTLINEARSOLVERPARAMETERS_HEADER_INCLUDED
#define OPM_PRINTLINEARSOLVERPARAMETERS_HEADER_INCLUDED

#include <ostream>
#include <vector>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm::detail
{

/// \brief Write a single property tree in JSON format to the stream
///
/// \note Helper function to write a property tree in JSON format to the log.
inline void writeJsonToStream(const Opm::PropertyTree& prm, std::ostream& os)
{
    // Write the property tree in JSON format to the output stream
    prm.write_json(os, true);
}

/// \brief Write a vector of property trees in JSON format to the stream.
///
/// \note Helper function to write a property tree in JSON format to the log.
inline void writeJsonToStream(const std::vector<Opm::PropertyTree>& prms, std::ostream& os)
{
    // Write each property tree in the vector in JSON format to the output stream
    for (const auto& p : prms) {
        writeJsonToStream(p, os);
    }
}

/// \brief Print the linear solver parameters to the log if requested.
///
/// \param parameters The linear solver parameters.
/// \param prm The property tree containing the parameters.
/// \param comm The communication object to check if we are on the IO rank.
///
/// This function will only print the parameters if the `linear_solver_print_json_definition_`
/// flag is set to true in parameters. It will print the prm in JSON format.
///
/// \tparam Comm The type of the communication object, which should support `rank()` method.
template <class VectorOrSingle, class Comm>
void
printLinearSolverParameters(const FlowLinearSolverParameters& parameters,
                            const VectorOrSingle& prm,
                            const Comm& comm)
{
    // Check if we are on the IO rank
    const bool on_io_rank = comm.rank() == 0;
    if (on_io_rank && parameters.linear_solver_print_json_definition_) {
        std::ostringstream os;
        os << "\nProperty tree for linear solvers:\n";
        writeJsonToStream(prm, os);
        // Write the parameters in JSON format
        OpmLog::note(os.str());
    }
}

/// \brief Print the linear solver parameters to the log if requested.
///
/// \param parameters The vector of linear solver parameters.
/// \param activeSolverNum The index of the active solver.
/// \param prm The vector property tree containing the parameters.
/// \param comm The communication object to check if we are on the IO rank.
///
/// This function will only print the parameters if the `linear_solver_print_json_definition_`
/// flag is set to true for the active solver. It will print the prm in JSON format.
///
/// \tparam Comm The type of the communication object, which should support `rank()` method.
template <class Comm>
void
printLinearSolverParameters(const std::vector<FlowLinearSolverParameters>& parameters,
                            int activeSolverNum,
                            const std::vector<Opm::PropertyTree>& prm,
                            const Comm& comm)
{
    printLinearSolverParameters(parameters[activeSolverNum], prm, comm);
}

} // namespace Opm::detail

#endif // OPM_PRINTLINEARSOLVERPARAMETERS_HEADER_INCLUDED
