#include <cuda.h>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/cusparse_safe_call.hpp>

namespace Opm::cuistl
{

template <class T>
CuSparseMatrix<T>::CuSparseMatrix(const T* nonZeroElements,
                                  const int* rowIndices,
                                  const int* columnIndices,
                                  int numberOfNonzeroElements,
                                  int blockSize,
                                  int numberOfRows)
    : nonZeroElements(nonZeroElements, numberOfNonzeroElements)
    , rowIndices(rowIndices, numberOfRows + 1)
    , columnIndices(columnIndices, numberOfNonzeroElements / blockSize)
    , numberOfNonzeroElements(numberOfNonzeroElements)
    , numberOfRows(numberOfRows)
{
    OPM_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&matrixDescription));
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(matrixDescription, baseType));

    // TODO: Are there better options for our matrices? Probably not.
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatType(matrixDescription, CUSPARSE_MATRIX_TYPE_GENERAL));
}

template <class T>
CuSparseMatrix<T>::~CuSparseMatrix()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(matrixDescription));
}

template <typename T>
void
CuSparseMatrix<T>::setUpperTriangular()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(matrixDescription, CUSPARSE_FILL_MODE_UPPER));
}

template <typename T>
void
CuSparseMatrix<T>::setLowerTriangular()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(matrixDescription, CUSPARSE_FILL_MODE_LOWER));
}

template <typename T>
void
CuSparseMatrix<T>::setUnitDiagonal()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(matrixDescription, CUSPARSE_DIAG_TYPE_UNIT));
}

template <typename T>
void
CuSparseMatrix<T>::setNonUnitDiagonal()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(matrixDescription, CUSPARSE_DIAG_TYPE_NON_UNIT));
}
} // namespace Opm::cuistl
