/*
  Copyright 2020 Equinor ASA

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

#ifndef FPGA_BILU0_HEADER_INCLUDED
#define FPGA_BILU0_HEADER_INCLUDED

#include <vector>

#include <opm/simulators/linalg/bda/ILUReorder.hpp>
#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

namespace Opm
{
namespace Accelerator
{

/*
 * This class implements a Blocked ILU0 preconditioner, with output data
 * specifically formatted for the FPGA.
 * The decomposition and reorders of the rows of the matrix are done on CPU.
 */

template <unsigned int block_size>
class FPGABILU0
{

private:
    int N;       // number of rows of the matrix
    int Nb;      // number of blockrows of the matrix
    int nnz;     // number of nonzeroes of the matrix (scalar)
    int nnzbs;   // number of blocks of the matrix
    std::unique_ptr<BlockedMatrix<block_size> > LMat = nullptr, UMat = nullptr, LUMat = nullptr;
    std::shared_ptr<BlockedMatrix<block_size> > rMat = nullptr; // reordered mat
    double *invDiagVals = nullptr;
    std::vector<int> diagIndex;
    std::vector<int> toOrder, fromOrder;
    std::vector<int> rowsPerColor;
    int numColors;
    int verbosity;

    // sizes and arrays used during RDF generation
    std::vector<std::vector<double> > nnzValues, LnnzValues, UnnzValues;
    std::vector<short int> colIndices, LColIndices, UColIndices;
    std::vector<unsigned char> NROffsets, LNROffsets, UNROffsets;
    std::vector<int> PIndicesAddr, LPIndicesAddr, UPIndicesAddr;
    std::vector<int> colorSizes, LColorSizes, UColorSizes;
    std::vector<int> nnzValsSizes, LnnzValsSizes, UnnzValsSizes;
    std::vector<std::vector<int> > colIndicesInColor, LColIndicesInColor, UColIndicesInColor;

    int rowSize, valSize;
    int LRowSize, LValSize, LNumColors;
    int URowSize, UValSize, UNumColors;
    std::vector<double> blockDiag;
    ILUReorder opencl_ilu_reorder;
    bool level_scheduling = false, graph_coloring = false;
    int numResultPointers = 21;
    std::vector<void *> resultPointers;
    int numResultSizes = 18;
    std::vector<int> resultSizes;
    int maxRowsPerColor, maxColsPerColor, maxNNZsPerRow, maxNumColors; // are set via the constructor

public:

    FPGABILU0(ILUReorder opencl_ilu_reorder, int verbosity, int maxRowsPerColor, int maxColsPerColor, int maxNNZsPerRow, int maxNumColors);

    ~FPGABILU0();

    // analysis (optional)
    bool init(BlockedMatrix<block_size> *mat);

    // ilu_decomposition
    bool create_preconditioner(BlockedMatrix<block_size> *mat);

    int* getToOrder()
    {
        return toOrder.data();
    }

    int* getFromOrder()
    {
        return fromOrder.data();
    }

    BlockedMatrix<block_size>* getRMat()
    {
        return rMat.get();
    }

    void **getResultPointers()
    {
        return resultPointers.data();
    }

    int *getResultSizes()
    {
        return resultSizes.data();
    }

};

} // namespace Accelerator
} // namespace Opm

#endif // FPGA_BILU0_HEADER_INCLUDED
