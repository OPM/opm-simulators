/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2019 Dr. Blatt - HPC-Simulation-Software & Services.

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

#ifndef OPM_MATRIXMARKETUTILS_HEADER_INCLUDED
#define OPM_MATRIXMARKETUTILS_HEADER_INCLUDED

#include <dune/istl/matrixmarket.hh>
#include <cstddef>

namespace Dune
{

namespace MatrixMarketImpl
{
    template <typename T, typename A, int brows, int bcols, typename D>
    void makeSparseEntries(Dune::BCRSMatrix<Dune::FieldMatrix<T, brows, bcols>, A>& matrix,
                           std::vector<int>& i,
                           std::vector<int>& j,
                           std::vector<double>& val,
                           const D&)
    {
        // addapted from readSparseEntries
        typedef Dune::BCRSMatrix<Dune::FieldMatrix<T, brows, bcols>, A> Matrix;
        std::vector<std::set<IndexData<D>>> rows(matrix.N() * brows);
        for (std::size_t kk = 0; kk < i.size(); ++kk) {
            std::size_t row;
            IndexData<D> data;
            row = i[kk];
            assert(row / bcols < matrix.N());
            data.number = val[kk];
            data.index = j[kk];
            assert(data.index / bcols < matrix.M());
            rows[row].insert(data);
        }

        // Setup the matrix sparsity pattern
        int nnz = 0;
        for (typename Matrix::CreateIterator iter = matrix.createbegin(); iter != matrix.createend(); ++iter) {
            for (std::size_t brow = iter.index() * brows, browend = iter.index() * brows + brows; brow < browend;
                 ++brow) {
                typedef typename std::set<IndexData<D>>::const_iterator Siter;
                for (Siter siter = rows[brow].begin(), send = rows[brow].end(); siter != send; ++siter, ++nnz)
                    iter.insert(siter->index / bcols);
            }
        }

        // Set the matrix values
        matrix = 0;

        MatrixValuesSetter<D, brows, bcols> Setter;

        Setter(rows, matrix);
    }
} // end namespace MatrixMarketImpl


template <typename T, typename A, int brows, int bcols>
void
makeMatrixMarket(Dune::BCRSMatrix<Dune::FieldMatrix<T, brows, bcols>, A>& matrix,
                 std::vector<int> i,
                 std::vector<int> j,
                 std::vector<T> val,
                 size_t rows,
                 size_t cols)
{
    // addapted from readMatrixMarket
    using namespace MatrixMarketImpl;
    // std::size_t rows, cols, entries;
    // std::size_t nnz, blockrows, blockcols;
    // std::tie(blockrows, blockcols, nnz) = calculateNNZ<brows, bcols>(rows, cols, entries, header);
    std::size_t blockrows = rows / brows;
    std::size_t blockcols = cols / bcols;
    matrix.setSize(blockrows, blockcols);
    matrix.setBuildMode(Dune::BCRSMatrix<Dune::FieldMatrix<T, brows, bcols>, A>::row_wise);
    makeSparseEntries(matrix, i, j, val, NumericWrapper<T>());
}


} // end namespace Dune



#endif // OPM_MATRIXMARKETUTILS_HEADER_INCLUDED
