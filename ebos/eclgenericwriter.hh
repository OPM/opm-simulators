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
 *
 * \copydoc Opm::EclWriter
 */
#ifndef EWOMS_ECL_GENERIC_WRITER_HH
#define EWOMS_ECL_GENERIC_WRITER_HH

#include "collecttoiorank.hh"
#include <ebos/ecltransmissibility.hh>

#include <opm/models/parallel/tasklets.hh>

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

namespace Action { class State; }
class EclipseIO;
class EclipseState;
class Inplace;
struct NNCdata;
class Schedule;
class SummaryConfig;
class SummaryState;
class UDQState;
class WellTestState;

template <class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
class EclGenericWriter
{
    using CollectDataToIORankType = CollectDataToIORank<Grid,EquilGrid,GridView>;
    using TransmissibilityType = EclTransmissibility<Grid,GridView,ElementMapper,Scalar>;

public:
    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    EclGenericWriter(const Schedule& schedule,
                     const EclipseState& eclState,
                     const SummaryConfig& summaryConfig,
                     const Grid& grid,
                     const EquilGrid* equilGrid,
                     const GridView& gridView,
                     const Dune::CartesianIndexMapper<Grid>& cartMapper,
                     const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper,
                     bool enableAsyncOutput,
                     bool enableEsmry);

    const EclipseIO& eclIO() const;

    void writeInit();

    void setTransmissibilities(const TransmissibilityType* globalTrans)
    {
        globalTrans_ = globalTrans;
    }

protected:
    const TransmissibilityType& globalTrans() const;

    void doWriteOutput(const int                     reportStepNum,
                       const bool                    isSubStep,
                       data::Solution&&              localCellData,
                       data::Wells&&                 localWellData,
                       data::GroupAndNetworkValues&& localGroupAndNetworkData,
                       data::Aquifers&&              localAquiferData,
                       const Action::State& actionState,
                       const WellTestState& wtestState,
                       const UDQState& udqState,
                       const SummaryState& summaryState,
                       const std::vector<Scalar>& thresholdPressure,
                       Scalar curTime,
                       Scalar nextStepSize,
                       bool doublePrecision);

    void evalSummary(int reportStepNum,
                     Scalar curTime,
                     const std::map<std::size_t, double>& wbpData,
                     const data::Wells& localWellData,
                     const data::GroupAndNetworkValues& localGroupAndNetworkData,
                     const std::map<int,data::AquiferData>& localAquiferData,
                     const std::map<std::pair<std::string, int>, double>& blockData,
                     const std::map<std::string, double>& miscSummaryData,
                     const std::map<std::string, std::vector<double>>& regionData,
                     SummaryState& summaryState,
                     UDQState& udqState,
                     const Inplace& inplace,
                     const Inplace& initialInPlace);

    CollectDataToIORankType collectToIORank_;
    const Grid& grid_;
    const GridView& gridView_;
    const Schedule& schedule_;
    const EclipseState& eclState_;
    const SummaryConfig& summaryConfig_;
    std::unique_ptr<EclipseIO> eclIO_;
    std::unique_ptr<TaskletRunner> taskletRunner_;
    Scalar restartTimeStepSize_;
    const TransmissibilityType* globalTrans_ = nullptr;
    const Dune::CartesianIndexMapper<Grid>& cartMapper_;
    const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper_;
    const EquilGrid* equilGrid_;
    std::vector<std::size_t> wbp_index_list_;

private:
    data::Solution computeTrans_(const std::unordered_map<int,int>& cartesianToActive) const;
    std::vector<NNCdata> exportNncStructure_(const std::unordered_map<int,int>& cartesianToActive) const;
};

} // namespace Opm

#endif
