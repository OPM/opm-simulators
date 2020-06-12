/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_WRITESYSTEMMATRIXHELPER_HEADER_INCLUDED
#define OPM_WRITESYSTEMMATRIXHELPER_HEADER_INCLUDED

#include <dune/istl/matrixmarket.hh>
#include <opm/simulators/linalg/MatrixMarketSpecializations.hpp>


namespace Opm
{
namespace Helper
{
    template <class SimulatorType, class MatrixType, class VectorType, class Communicator>
    void writeSystem(const SimulatorType& simulator,
                     const MatrixType& matrix,
                     const VectorType& rhs,
                     [[maybe_unused]] const Communicator* comm,
                     std::string ename = std::string("")
        )
    {
        std::string dir = simulator.problem().outputDir();
        if (dir == ".") {
            dir = "";
        } else if (!dir.empty() && dir.back() != '/') {
            dir += "/";
        }
        namespace fs = Opm::filesystem;
        fs::path output_dir(dir);
        fs::path subdir("reports");
        output_dir = output_dir / subdir;
        if (!(fs::exists(output_dir))) {
            fs::create_directory(output_dir);
        }
        // Combine and return.
        std::ostringstream oss;
        oss << "prob_" << simulator.episodeIndex() << "_time_";
        oss << std::setprecision(15) << std::setw(12) << std::setfill('0') << simulator.time() << "_";
        int nit = simulator.model().newtonMethod().numIterations();
        oss << "_nit_" << nit << "_";
        std::string output_file(oss.str());
        fs::path full_path = output_dir / output_file;
        std::string prefix = full_path.string();
        {
            std::string filename = prefix + ename + "matrix_istl";
#if HAVE_MPI
            if (comm != nullptr) { // comm is not set in serial runs
                Dune::storeMatrixMarket(matrix, filename, *comm, true);
            } else
#endif
            {
                Dune::storeMatrixMarket(matrix, filename + ".mm");
            }
        }
        {
            std::string filename = prefix + ename + "rhs_istl";
#if HAVE_MPI
            if (comm != nullptr) { // comm is not set in serial runs
                Dune::storeMatrixMarket(rhs, filename, *comm, true);
            } else
#endif
            {
                Dune::storeMatrixMarket(rhs, filename + ".mm");
            }
        }
    }


} // namespace Helper
} // namespace Opm
#endif
