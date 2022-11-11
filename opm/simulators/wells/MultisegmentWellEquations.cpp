/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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

#include <config.h>
#include <opm/simulators/wells/MultisegmentWellEquations.hpp>

#include <opm/simulators/wells/MultisegmentWellGeneric.hpp>

namespace Opm {

template<class Scalar, int numWellEq, int numEq>
MultisegmentWellEquations<Scalar,numWellEq,numEq>::
MultisegmentWellEquations(const MultisegmentWellGeneric<Scalar>& well)
    : well_(well)
{
}

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::
init(const int num_cells,
     const int numPerfs,
     const std::vector<int>& cells)
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise);
    duneD_.setBuildMode(DiagMatWell::row_wise);

    // set the size and patterns for all the matrices and vectors
    // [A C^T    [x    =  [ res
    //  B D] x_well]      res_well]

    // calculating the NNZ for duneD_
    // NNZ = number_of_segments + 2 * (number_of_inlets / number_of_outlets)
    {
        int nnz_d = well_.numberOfSegments();
        for (const std::vector<int>& inlets : well_.segmentInlets()) {
            nnz_d += 2 * inlets.size();
        }
        duneD_.setSize(well_.numberOfSegments(), well_.numberOfSegments(), nnz_d);
    }
    duneB_.setSize(well_.numberOfSegments(), num_cells, numPerfs);
    duneC_.setSize(well_.numberOfSegments(), num_cells, numPerfs);

    // we need to add the off diagonal ones
    for (auto row = duneD_.createbegin(),
              end = duneD_.createend(); row != end; ++row) {
        // the number of the row corrspnds to the segment now
        const int seg = row.index();
        // adding the item related to outlet relation
        const Segment& segment = well_.segmentSet()[seg];
        const int outlet_segment_number = segment.outletSegment();
        if (outlet_segment_number > 0) { // if there is a outlet_segment
            const int outlet_segment_index = well_.segmentNumberToIndex(outlet_segment_number);
            row.insert(outlet_segment_index);
        }

        // Add nonzeros for diagonal
        row.insert(seg);

        // insert the item related to its inlets
        for (const int& inlet : well_.segmentInlets()[seg]) {
            row.insert(inlet);
        }
    }

    // make the C matrix
    for (auto row = duneC_.createbegin(),
              end = duneC_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : well_.segmentPerforations()[row.index()]) {
            const int cell_idx = cells[perf];
            row.insert(cell_idx);
        }
    }

    // make the B^T matrix
    for (auto row = duneB_.createbegin(),
              end = duneB_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : well_.segmentPerforations()[row.index()]) {
            const int cell_idx = cells[perf];
            row.insert(cell_idx);
        }
    }

    resWell_.resize(well_.numberOfSegments());
}

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::clear()
{
    duneB_ = 0.0;
    duneC_ = 0.0;
    duneD_ = 0.0;
    resWell_ = 0.0;
    duneDSolver_.reset();
}

#define INSTANCE(numWellEq, numEq) \
template class MultisegmentWellEquations<double,numWellEq,numEq>;

INSTANCE(2,1)
INSTANCE(2,2)
INSTANCE(2,6)
INSTANCE(3,2)
INSTANCE(3,3)
INSTANCE(3,4)
INSTANCE(4,3)
INSTANCE(4,4)
INSTANCE(4,5)

}
