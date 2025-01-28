/*
  Copyright 2024, SINTEF Digital

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

namespace Opm
{

template <typename Scalar, int numWellEq, int numEq>
CompWellEquations<Scalar, numWellEq, numEq>::
CompWellEquations()
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise),
    duneD_.setBuildMode(DiagMatWell::row_wise);
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}


template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
init(const int num_conn)
{
    duneD_.setSize(1, 1, 1);
    duneB_.setSize(1, num_conn, num_conn);
    duneC_.setSize(1, num_conn, num_conn);

    for (auto row = duneD_.createbegin(),
              end = duneD_.createend(); row != end; ++row) {
        // Add nonzeros for diagonal
        row.insert(row.index());
    }

    for (auto row = duneB_.createbegin(),
                 end = duneB_.createend(); row != end; ++row) {
        for (int con = 0 ; con < num_conn; ++con) {
            row.insert(con);
        }
    }

    // make the C^T matrix
    // TODO: let us see whether we should change the naming of DuneC_
    for (auto row = duneC_.createbegin(),
                 end = duneC_.createend(); row != end; ++row) {
        for (int con = 0; con < num_conn; ++con) {
            row.insert(con);
        }
    }

    resWell_.resize(1);

    // some others in the future
}

template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
clear()
{
    duneB_ = 0.0;
    duneC_ = 0.0;
    duneD_ = 0.0;
    resWell_ = 0.0;
}


} // end of namespace Opm