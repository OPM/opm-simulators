// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=4:
/*
  Copyright 2021 Equinor AS.

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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/**
 * \file
 * \brief A datahandle sending data located in multiple vectors
 * \author Markus Blatt, OPM-OP AS
 */

#include <dune/grid/common/datahandleif.hh>

#include <cstddef>

namespace Opm
{

/// \brief A data handle sending multiple data store in vectors attached
///        to cells.
///
/// Each data is assumed to a container with operator[] and the
/// class operates on a vector of these.
/// \tparam GridView the type of the grid view the data associated with
/// \tparam The type of the vector of vectors.
template<class GridView, class Vector>
class VectorVectorDataHandle
  : public Dune::CommDataHandleIF<VectorVectorDataHandle<GridView,Vector>,
                                  std::decay_t<decltype(Vector()[0][0])>>
{
public:

  /// \brief the data type we send
  using DataType = std::decay_t<decltype(Vector()[0][0])>;

  /// \brief Constructor
  /// \param data The vector of data vectors
  /// \param gridView The gridview the data is attached to.
  VectorVectorDataHandle(Vector& data, const GridView& gridView)
    : data_(data), gridView_(gridView)
  {}

  bool contains(int /* dim */, int codim) const
  {
    return codim == 0;
  }

#if DUNE_VERSION_LT(DUNE_GRID, 2, 8)
  bool fixedsize(int /* dim */, int /* codim */) const
  {
    return true;
  }
#else

  bool fixedSize(int /* dim */, int /* codim */) const
  {
    return true;
  }
#endif
  template<class EntityType>
  std::size_t size(const EntityType /* entity */) const
  {
    return data_.size();
  }


  template<class BufferType, class EntityType>
  void gather(BufferType& buffer, const EntityType& e) const
  {
    for(const auto& vec: data_)
    {
      buffer.write(vec[gridView_.indexSet().index(e)]);
    }
  }

  template<class BufferType, class EntityType>
  void scatter(BufferType& buffer, const EntityType& e,
               [[maybe_unused]] std::size_t n)
  {
    assert(n == data_.size());
    for(auto& vec: data_)
    {
      buffer.read(vec[gridView_.indexSet().index(e)]);
    }
  }
private:
  Vector& data_;
  const GridView& gridView_;
};

} // end namespace Opm
