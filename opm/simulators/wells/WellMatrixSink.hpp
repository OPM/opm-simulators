/*
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2017 Statoil ASA.

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

#ifndef OPM_WELLMATRIXSINK_HPP
#define OPM_WELLMATRIXSINK_HPP

#include <cstddef>
#include <functional>
#include <utility>

namespace Opm {

template<class Block>
class WellMatrixSink
{
public:
    using MatrixBlock = Block;
    using AddBlock = std::function<void(std::size_t, std::size_t, const Block&)>;

    explicit WellMatrixSink(AddBlock addBlock)
        : addBlock_(std::move(addBlock))
    {}

    void addToBlock(std::size_t row, std::size_t column, const Block& value) const
    {
        addBlock_(row, column, value);
    }

private:
    AddBlock addBlock_;
};

} // namespace Opm

#endif
