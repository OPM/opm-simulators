
/*
  Copyright 2022-2023 SINTEF AS
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

#ifndef OPM_DILUUTIL_HEADER_INCLUDED
#define OPM_DILUUTIL_HEADER_INCLUDED

#include <config.h>
#include <fstream>
#include <mpi.h>
#include <opm/grid/utility/SparseTable.hpp>

namespace Opm::DILUUtils{

// TODO: make proper doxygen
// This function is intended to be used by parallel DILU implementations to
// explore how parallelizable different linear systems are. The results
// are for now written to a file from where flow is run
template <class T>
void writeSparseTableRowSizesToFile(const Opm::SparseTable<T> *sparseTable){

        int rank = 0;
        int size = sparseTable->size();
#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        std::string filename = "DILU_parallelism_distribution_on_rank_" + std::to_string(rank) + ".txt";
        std::ofstream file(filename);

        file << "PRINTING NUMBER OF LEVEL SETS, THEN NUMBER OF LEVELS IN EACH SET" << std::endl;
        file << size << std::endl;

        for (int i = 0; i < size; i++){
            file << sparseTable->rowSize(i) << std::endl;
        }
}

template void writeSparseTableRowSizesToFile(const Opm::SparseTable<size_t>*);

} // END NAMESPACE OPM

#endif
