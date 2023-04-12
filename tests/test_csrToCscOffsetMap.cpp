#include <config.h>
#include <fstream>

#define BOOST_TEST_MODULE CsrToCscOffsetMap
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>

#include <opm/simulators/linalg/bda/opencl/BISAI.hpp>

BOOST_AUTO_TEST_CASE(testcsrtocscoffsetmap){

    using Matrix = Dune::BCRSMatrix<double>;

    Matrix matrix;
    {
        std::ifstream mfile("offset_map_matrix.txt");
        if (!mfile) {
            throw std::runtime_error("Could not read matrix file");
        }
        readMatrixMarket(matrix, mfile);
    }

    // a transposed version of the matrix is read because the transposed
    // of a CSR representation is equivalente to CSC
    Matrix matrix_transposed;
    {
        std::ifstream mfile("offset_map_matrix_transposed.txt");
        if (!mfile) {
            throw std::runtime_error("Could not read matrix file");
        }
        readMatrixMarket(matrix_transposed, mfile);
    }

    // has to make copy because the output of readMatrixMarket does not
    // have contiguous non-zero values
    Matrix matrix_copy(matrix);
    Matrix matrix_transposed_copy(matrix_transposed);

    std::vector<int> rowPointers, colIndices, map;

    auto* nnzValues = &matrix_copy[0][0];
    auto* nnzValues_transposed = &matrix_transposed_copy[0][0];

    rowPointers.emplace_back(0);
    for (Matrix::Iterator r = matrix.begin(); r != matrix.end(); ++r) {
        for (auto c = r->begin(); c != r->end(); ++c) {
            colIndices.emplace_back(c.index());
        }
        rowPointers.emplace_back(colIndices.size());
    }

    map = Opm::Accelerator::buildCsrToCscOffsetMap(rowPointers, colIndices);

    for (unsigned int i = 0; i < colIndices.size(); i++){
        BOOST_CHECK_EQUAL(nnzValues[i], nnzValues_transposed[map[i]]);
    }
}
