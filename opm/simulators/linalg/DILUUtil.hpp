
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
/// @brief This function writes to a file a sparse table and the sizes of each row
/// The intended use is to be called from a DILU preconditioner to explore how
/// parallelizable the upper and lower solves are. For now the function creats a file
/// in the directory from which flow was called, only intended to be done when verbosity>0
/// @tparam T Any type, it does not matter what objects are stored in the table
/// @param sparse_table A pointer to a SparseTable, this function does not change the object
template <class T>
void writeSparseTableRowSizesToFile(const Opm::SparseTable<T> *sparse_table){

        int rank = 0;
        int size = sparse_table->size();
#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        std::string filename = "DILU_parallelism_distribution_on_rank_" + std::to_string(rank) + ".txt";
        std::ofstream file(filename);

        // print brief exeplenation of how to interpret the data in the file
        file << "First is the number of levels, then the size of each level set" << std::endl;
        file << size << std::endl;

        for (int i = 0; i < size; i++){
            file << sparse_table->rowSize(i) << std::endl;
        }
}

template void writeSparseTableRowSizesToFile(const Opm::SparseTable<size_t>*);

} // namespace Opm::DILUUtils

#endif
