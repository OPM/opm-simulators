/*
  Copyright SINTEF AS

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
#ifndef OPM_CUVECTOR_HEADER_INCLUDED
#define OPM_CUVECTOR_HEADER_INCLUDED
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <exception>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuBlasHandle.hpp>
#include <vector>

namespace Opm::cuistl
{

//! \brief Simple vector class on the GPU.
template <typename T>
class CuVector
{
public:
    using field_type = T;
    using size_type = size_t;

    CuVector(const CuVector<T>& other);
    explicit CuVector(const std::vector<T>& data);
    CuVector& operator=(const CuVector<T>&);
    CuVector& operator=(T scalar);

    explicit CuVector(const int numberOfElements);
    CuVector(const T* dataOnHost, const int numberOfElements);
    virtual ~CuVector();

    T* data();
    const T* data() const;

    template <int BlockDimension>
    void copyFromHost(const Dune::BlockVector<Dune::FieldVector<T, BlockDimension>>& vector)
    {
        if (m_numberOfElements != vector.dim()) {
            OPM_THROW(std::runtime_error,
                      "Given incompatible vector size. CuVector has size " + std::to_string(m_numberOfElements)
                          + ",\nhowever, BlockVector has N() = " + std::to_string(vector.N())
                          + ", and dim() = " + std::to_string(vector.dim()));
        }
        const auto dataPointer = static_cast<const T*>(&(vector[0][0]));
        copyFromHost(dataPointer, m_numberOfElements);
    }

    template <int BlockDimension>
    void copyToHost(Dune::BlockVector<Dune::FieldVector<T, BlockDimension>>& vector) const
    {
        if (m_numberOfElements != vector.dim()) {
            OPM_THROW(std::runtime_error,
                      "Given incompatible vector size. CuVector has size " + std::to_string(m_numberOfElements)
                          + ",\nhowever, the BlockVector has has N() = " + std::to_string(vector.N())
                          + ", and dim() = " + std::to_string(vector.dim()));
        }
        const auto dataPointer = static_cast<T*>(&(vector[0][0]));
        copyToHost(dataPointer, m_numberOfElements);
    }
    void copyFromHost(const T* dataPointer, int numberOfElements);
    void copyToHost(T* dataPointer, int numberOfElements) const;
    void copyFromHost(const std::vector<T>& data);
    void copyToHost(std::vector<T>& data) const;

    CuVector<T>& operator*=(const T& scalar);

    CuVector<T>& axpy(T alpha, const CuVector<T>& y);
    CuVector<T>& operator+=(const CuVector<T>& other);
    CuVector<T>& operator-=(const CuVector<T>& other);
    T dot(const CuVector<T>& other) const;
    T two_norm() const;

    //! Computes the dot product sum_i this[indexSet[i]] * other[indexSet[i]]
    T dot(const CuVector<T>& other, const CuVector<int>& indexSet, CuVector<T>& buffer) const;
    //! Computes the norm sqrt(sum_i this[indexSet[i]] * this[indexSet[i]])
    T two_norm(const CuVector<int>& indexSet, CuVector<T>& buffer) const;
    //! Computes the dot product sum_i this[indexSet[i]] * other[indexSet[i]]
    T dot(const CuVector<T>& other, const CuVector<int>& indexSet) const;
    //! Computes the norm sqrt(sum_i this[indexSet[i]] * this[indexSet[i]])
    T two_norm(const CuVector<int>& indexSet) const;


    size_type dim() const;

    std::vector<T> asStdVector() const;

    template <int blockSize>
    Dune::BlockVector<Dune::FieldVector<T, blockSize>> asDuneBlockVector() const
    {
        OPM_ERROR_IF(dim() % blockSize != 0,
                     "blockSize is not a multiple of dim(). Given blockSize = " << blockSize
                                                                                << " and dim() = " << dim());

        Dune::BlockVector<Dune::FieldVector<T, blockSize>> returnValue(dim() / blockSize);
        copyToHost(returnValue);
        return returnValue;
    }


    void setZeroAtIndexSet(const CuVector<int>& indexSet);

private:
    T* m_dataOnDevice = nullptr;
    const int m_numberOfElements;
    impl::CuBlasHandle& m_cuBlasHandle;
};
} // namespace Opm::cuistl
#endif
