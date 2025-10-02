//===========================================================================
//
// File: AluGridLevelCartesianIndexMapper.hpp
//
// Created: Tue October 01  09:44:00 2024
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2024 Equinor ASA.

  This file is part of The Open Porous Media project  (OPM).

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
#ifndef OPM_ALUGRIDLEVELCARTESIANINDEXMAPPER_HPP
#define OPM_ALUGRIDLEVELCARTESIANINDEXMAPPER_HPP

#include <dune/alugrid/grid.hh>
#include <opm/grid/common/LevelCartesianIndexMapper.hpp>
#include <opm/simulators/flow/AluGridCartesianIndexMapper.hpp>

#include <array>
#include <memory>

namespace Opm {

// Interface class to access the local Cartesian grid of each level grid (when refinement).
// Further documentation in opm/grid/common/LevelCartesianIndexMapper.hpp
//
// Adapter Design Pattern: In this case, LevelCartesianIndexMapper uses the Object Adapter variant, where it holds an instance
// (here, a std::unique_ptr) of CartesianIndexMapper, the wrapped type. The goal is to provide a standardized interface, allowing
// incompatible functionality (such as Cartesian indexing in the context of refinement that may not be supported - yet -for all
// grid types, like CpGrid) to integrate smoothly within the existing conventions.
//
// Specialization for AluGrid
template<>
class LevelCartesianIndexMapper<Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>>
{

#if HAVE_MPI
    using Grid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using Grid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI

 public:
    static constexpr int dimension = 3 ;

    explicit LevelCartesianIndexMapper(const Dune::CartesianIndexMapper<Grid>& cartesianIndexMapper)
        : cartesianIndexMapper_{std::make_unique<Dune::CartesianIndexMapper<Grid>>(cartesianIndexMapper)}
    {}

    const std::array<int,3>& cartesianDimensions(int level) const
    {
        throwIfLevelPositive(level);
        return cartesianIndexMapper_ ->cartesianDimensions();
    }

    int cartesianSize(int level) const
    {
        throwIfLevelPositive(level);
        return cartesianIndexMapper_->cartesianSize();
    }

    int compressedSize(int level) const
    {
        throwIfLevelPositive(level);
        return cartesianIndexMapper_-> compressedSize();
    }

    int cartesianIndex( const int compressedElementIndex, const int level) const
    {
        throwIfLevelPositive(level);;
        return cartesianIndexMapper_->cartesianIndex(compressedElementIndex);
    }

    void cartesianCoordinate(const int compressedElementIndex, std::array<int,dimension>& coords, int level) const
    {
        throwIfLevelPositive(level);
        cartesianIndexMapper_->cartesianCoordinate(compressedElementIndex, coords);
    }

 private:
    std::unique_ptr<Dune::CartesianIndexMapper<Grid>> cartesianIndexMapper_;

    void throwIfLevelPositive(int level) const
    {
        if (level) {
            throw std::invalid_argument("Invalid level.\n");
        }
    }
};

}

#endif
