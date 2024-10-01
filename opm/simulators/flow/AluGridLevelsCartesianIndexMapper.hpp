//===========================================================================
//
// File: LevelsCartesianIndexMapper.hh
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
#ifndef OPM_POLYHEDRALGRIDLEVELSCARTESIANINDEXMAPPER_HH
#define OPM_POLYHEDRALGRIDLEVELSCARTESIANINDEXMAPPER_HH

#include <opm/grid/common/LevelsCartesianIndexMapper.hh>
#include <opm/grid/polyhedralgrid.hh>

#include <memory>


namespace Dune
{
template<...>
class AluGrid;
}

namespace Opm
{

template<...>
class LevelsCartesianIndexMapper<<Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>>
{
   
#if HAVE_MPI
    using Grid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>; 
#else    
    using Grid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI
    
public:
    static const int dimension = 3 ;

    explicit LevelsCartesianIndexMapper(const Grid& grid) : grid_{ &grid }// std::make_unique<Grid>(grid)}
    {}

    const std::array<int,3>& levelCartesianDimensions(int level) const
    {
        throwIfLevelPositive(level);
        return grid_->logicalCartesianSize();
    }

    int levelCartesianSize(int level) const
    {
        throwIfLevelPositive(level);
        return computeLevelCartesianSize(0);
    }

    int levelCompressedSize(int level) const
    {
        throwIfLevelPositive(level);
        return grid_->size(0);
    }


    int cartesianIndexLevel( const int compressedElementIndex, const int level) const
    {
        throwIfLevelPositive(level);
        assert( compressedElementIndex >= 0 && compressedElementIndex < levelCompressedSize(0) );
        return grid_->globalCell()[ compressedElementIndex ];
    }

    void cartesianCoordinateLevel(const int compressedElementIndex, std::array<int,dimension>& coords, int level) const
    {
        throwIfLevelPositive(level);

        int gc = cartesianIndexLevel( compressedElementIndex, 0);
        auto cartesianDimensions = grid_->logicalCartesianSize();
        if( dimension >=2 )
        {
            for( int d=0; d<dimension-2; ++d )
            {
                coords[d] = gc % cartesianDimensions[d];  gc /= cartesianDimensions[d];
            }

            coords[dimension-2] = gc % cartesianDimensions[dimension-2];
            coords[dimension-1] = gc / cartesianDimensions[dimension-1];
        }
        else
            coords[ 0 ] = gc ;
    }

private:
    const Grid* grid_;

    int computeLevelCartesianSize(int level) const
    {
        int size = levelCartesianDimensions(level)[ 0 ];
        for( int d=1; d<dimension; ++d )
            size *= levelCartesianDimensions(level)[ d ];
        return size;
    }

    void throwIfLevelPositive(int level) const
    {
        if (level) {
            throw std::invalid_argument("Invalid level.\n");
        }
    }
};

}

#endif
