/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILMODELENUMS_HEADER_INCLUDED
#define OPM_BLACKOILMODELENUMS_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>

#include <array>
#include <vector>

namespace Opm
{

    constexpr const auto Water        = BlackoilPhases::Aqua;
    constexpr const auto Oil          = BlackoilPhases::Liquid;
    constexpr const auto Gas          = BlackoilPhases::Vapour;
    constexpr const auto MaxNumPhases = BlackoilPhases::MaxNumPhases;

    enum PrimalVariables {
        Sg = 0,
        RS = 1,
        RV = 2
    };

    enum CanonicalVariablePositions {
        Pressure = 0,
        Sw = 1,
        Xvar = 2,
        Qs = 3,
        Bhp = 4,
        Next // For extension.
    };

    struct FIPDataEnums {

        enum FipId {
            FIP_AQUA = Opm::Water,
            FIP_LIQUID = Opm::Oil,
            FIP_VAPOUR = Opm::Gas,
            FIP_DISSOLVED_GAS = 3,
            FIP_VAPORIZED_OIL = 4,
            FIP_PV = 5,                    //< Pore volume
            FIP_WEIGHTED_PRESSURE = 6
        };

        static const int fipValues = FIP_WEIGHTED_PRESSURE + 1 ;
    };

    class FIPData : public FIPDataEnums
    {
    public:
        typedef std::vector<double> VectorType;

        using FIPDataEnums :: FipId;
        using FIPDataEnums :: fipValues ;
        std::array< VectorType, fipValues> fip;

        // default constructor
        FIPData() {}

        // initialize from array of Eigen vectors (or std::vectors)
        template <class V>
        explicit FIPData( const std::array< V, fipValues>& otherFip )
        {
            // copy fip vector from V to std::vector
            for( int i=0; i<fipValues; ++i ) {
                fip[ i ] = VectorType(otherFip[ i ].data(), otherFip[ i ].data() + otherFip[ i ].size() );
            }
        }
    };

} // namespace Opm

#endif // OPM_BLACKOILMODELENUMS_HEADER_INCLUDED
