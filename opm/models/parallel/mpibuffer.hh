// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::MpiBuffer
 */
#ifndef EWOMS_MPI_BUFFER_HH
#define EWOMS_MPI_BUFFER_HH

#if HAVE_MPI
#include <mpi.h>
#include <type_traits>
#endif

#include <cassert>
#include <cstddef>
#include <vector>

namespace Opm {

/*!
 * \brief Simplifies handling of buffers to be used in conjunction with MPI
 */
template <class DataType>
class MpiBuffer
{
public:
    MpiBuffer()
    {
        setMpiDataType_();
        updateMpiDataSize_();
    }

    explicit MpiBuffer(std::size_t size)
    {
        data_.resize(size);

        setMpiDataType_();
        updateMpiDataSize_();
    }

    /*!
     * \brief Set the size of the buffer
     */
    void resize(std::size_t newSize)
    {
        data_.resize(newSize);
        updateMpiDataSize_();
    }

    /*!
     * \brief Send the buffer asyncronously to a peer process.
     */
    void send([[maybe_unused]] unsigned peerRank)
    {
#if HAVE_MPI
        MPI_Isend(data_.data(),
                  static_cast<int>(mpiDataSize_),
                  mpiDataType_,
                  static_cast<int>(peerRank),
                  0, // tag
                  MPI_COMM_WORLD,
                  &mpiRequest_);
#endif
    }

    /*!
     * \brief Wait until the buffer was send to the peer completely.
     */
    void wait()
    {
#if HAVE_MPI
        MPI_Wait(&mpiRequest_, &mpiStatus_);
#endif // HAVE_MPI
    }

    /*!
     * \brief Receive the buffer syncronously from a peer rank
     */
    void receive([[maybe_unused]] unsigned peerRank)
    {
#if HAVE_MPI
        MPI_Recv(data_.data(),
                 static_cast<int>(mpiDataSize_),
                 mpiDataType_,
                 static_cast<int>(peerRank),
                 0, // tag
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
#endif // HAVE_MPI
    }

#if HAVE_MPI
    /*!
     * \brief Returns the current MPI_Request object.
     *
     * This object is only well defined after the send() method.
     */
    MPI_Request& request()
    { return mpiRequest_; }

    /*!
     * \brief Returns the current MPI_Request object.
     *
     * This object is only well defined after the send() method.
     */
    const MPI_Request& request() const
    { return mpiRequest_; }

    /*!
     * \brief Returns the current MPI_Status object.
     *
     * This object is only well defined after the receive() and wait() methods.
     */
    MPI_Status& status()
    { return mpiStatus_; }

    /*!
     * \brief Returns the current MPI_Status object.
     *
     * This object is only well defined after the receive() and wait() methods.
     */
    const MPI_Status& status() const
    { return mpiStatus_; }
#endif // HAVE_MPI

    /*!
     * \brief Returns the number of data objects in the buffer
     */
    std::size_t size() const
    { return data_.size(); }

    /*!
     * \brief Provide access to the buffer data.
     */
    DataType& operator[](std::size_t i)
    {
        assert(i < data_.size());
        return data_[i];
    }

    /*!
     * \brief Provide access to the buffer data.
     */
    const DataType& operator[](std::size_t i) const
    {
        assert(i < data_.size());
        return data_[i];
    }

private:
    void setMpiDataType_()
    {
#if HAVE_MPI
        // set the MPI data type
        if constexpr (std::is_same_v<DataType, char>) {
            mpiDataType_ = MPI_CHAR;
        }
        else if constexpr (std::is_same_v<DataType, unsigned char>) {
            mpiDataType_ = MPI_UNSIGNED_CHAR;
        }
        else if constexpr (std::is_same_v<DataType, short>) {
            mpiDataType_ = MPI_SHORT;
        }
        else if constexpr (std::is_same_v<DataType, unsigned short>) {
            mpiDataType_ = MPI_UNSIGNED_SHORT;
        }
        else if constexpr (std::is_same_v<DataType, int>) {
            mpiDataType_ = MPI_INT;
        }
        else if constexpr (std::is_same_v<DataType, unsigned>) {
            mpiDataType_ = MPI_UNSIGNED;
        }
        else if constexpr (std::is_same_v<DataType, long>) {
            mpiDataType_ = MPI_LONG;
        }
        else if constexpr (std::is_same_v<DataType, unsigned long>) {
            mpiDataType_ = MPI_UNSIGNED_LONG;
        }
        else if constexpr (std::is_same_v<DataType, long long>) {
            mpiDataType_ = MPI_LONG_LONG;
        }
        else if constexpr (std::is_same_v<DataType, unsigned long long>) {
            mpiDataType_ = MPI_UNSIGNED_LONG_LONG;
        }
        else if constexpr (std::is_same_v<DataType, float>) {
            mpiDataType_ = MPI_FLOAT;
        }
        else if constexpr (std::is_same_v<DataType, double>) {
            mpiDataType_ = MPI_DOUBLE;
        }
        else if constexpr (std::is_same_v<DataType, long double>) {
            mpiDataType_ = MPI_LONG_DOUBLE;
        }
        else {
            mpiDataType_ = MPI_BYTE;
        }
#endif // HAVE_MPI
    }

    void updateMpiDataSize_()
    {
#if HAVE_MPI
        mpiDataSize_ = data_.size();
        if (mpiDataType_ == MPI_BYTE) {
            mpiDataSize_ *= sizeof(DataType);
        }
#endif // HAVE_MPI
    }

    std::vector<DataType> data_;

#if HAVE_MPI
    std::size_t mpiDataSize_;
    MPI_Datatype mpiDataType_;
    MPI_Request mpiRequest_;
    MPI_Status mpiStatus_;
#endif // HAVE_MPI
};

} // namespace Opm

#endif
