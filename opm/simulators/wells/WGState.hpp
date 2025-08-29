/*
  Copyright 2021 Equinor ASA

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

#ifndef OPM_WGSTATE_HEADER_INCLUDED
#define OPM_WGSTATE_HEADER_INCLUDED

#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <memory>

namespace Opm {

template<class Scalar> class ParallelWellInfo;

/*
  Microscopic class to handle well, group and well test state.
*/

template<typename Scalar, typename IndexTraits>
struct WGState
{
    explicit WGState(const PhaseUsageInfo<IndexTraits>& pu);

    static WGState serializationTestObject(const ParallelWellInfo<Scalar>& pinfo);

    void wtest_state(std::unique_ptr<WellTestState> wtest_state);

    WellState<Scalar, IndexTraits> well_state;
    GroupState<Scalar> group_state;
    WellTestState well_test_state;

    bool operator==(const WGState&) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(well_state);
        serializer(group_state);
        serializer(well_test_state);
    }
};

}

#endif
