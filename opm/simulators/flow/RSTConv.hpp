/*
  Copyright 2023 SINTEF Digital

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

#ifndef OPM_RST_CONV_HEADER_INCLUDED
#define OPM_RST_CONV_HEADER_INCLUDED

#include <array>
#include <vector>

#include <opm/simulators/utils/ParallelCommunication.hpp>

namespace Opm {

class RSTConfig;

//! \brief Class computing RPTRST CONV output.
class RSTConv
{
public:
    //! \brief Constructor.
    //! \param globalCell Mapping from local to global cell indices
    //! \param comm Parallel communicator
    RSTConv(const std::vector<int>& globalCell,
            Parallel::Communication comm)
        : globalCell_(globalCell)
        , comm_(comm)
    {}

    //! \brief Init state at beginning of step.
    //! \param numCells Global number of active cells in the model
    //! \param rst_config RPTRST configuration
    //! \param compIdx Component index for phases {OIL, GAS, WAT}, -1 if inactive
    void init(const std::size_t numCells,
              const RSTConfig& rst_config,
              const std::array<int,3>& compIdx);

    //! \brief Adds the CONV output for given residual vector.
    template<class ResidualVector>
    void update(const ResidualVector& resid);

    //! \brief Obtain a const-ref to the accumulated data.
    const std::vector<std::vector<int>>& getData() const
    { return cnv_X_; }

private:
    //! \brief Gathers and accumulates results to the CONV arrays.
    //! \param lIdx Vector of local indices (N first is used)
    //! \param resid Residual vector
    //! \param comp Component to consider
    //! \details Handles parallel reduction
    template<class ResidualVector>
    void gatherAndAccumulate(const std::vector<int>& lIdx,
                             const ResidualVector& resid, int comp);

    const std::vector<int>& globalCell_; //!< Global cell indices
    Parallel::Communication comm_; //!< Communicator
    std::vector<std::vector<int>> cnv_X_; //!< Counts of worst cells for RPTRST CONV
    std::array<int,3> compIdx_; //!< Component indices
    int N_ = 0; //!< Number of cells to consider
};

} // namespace Opm

#endif // OPM_RST_CONV_HEADER_INCLUDED
