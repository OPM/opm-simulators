/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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

#ifndef OPM_BLACKOILWELLMODEL_RESTART_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_RESTART_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/output/data/Wells.hpp>

#include <vector>

namespace Opm {

class BlackoilWellModelGeneric;
namespace data {
struct GroupData;
class GroupAndNetworkValues;
}
class GroupState;
class GuideRate;
class GuideRateConfig;
struct PerforationData;
struct PhaseUsage;
class SingleWellState;
class WellState;

/// Class for restarting the blackoil well model.
class BlackoilWellModelRestart
{
public:
    //! \brief Constructor initializes reference to the well model.
    BlackoilWellModelRestart(const BlackoilWellModelGeneric& wellModel)
        : wellModel_(wellModel)
    {}

    //! \brief Loads guide rates from restart structures.
    void loadRestartGuideRates(const int                    report_step,
                               const GuideRateModel::Target target,
                               const data::Wells&           rst_wells,
                               GuideRate&                   guide_rate) const;

    //! \brief Loads guide rates from restart structures.
    void loadRestartGuideRates(const int                                     report_step,
                               const GuideRateConfig&                        config,
                               const std::map<std::string, data::GroupData>& rst_groups,
                               GuideRate&                                    guide_rate) const;

    //! \brief Loads well data from restart structures.
    void loadRestartData(const data::Wells&                 rst_wells,
                         const data::GroupAndNetworkValues& grpNwrkValues,
                         const bool                         handle_ms_well,
                         WellState&                         well_state,
                         GroupState&                        grpState) const;

private:
    //! \brief Loads per-connection data from restart structures.
    void loadRestartConnectionData(const std::vector<data::Rates::opt>& phs,
                                   const data::Well&                    rst_well,
                                   const std::vector<PerforationData>&  old_perf_data,
                                   SingleWellState&                     ws) const;

    //! \brief Loads per-segment data from restart structures.
    void loadRestartSegmentData(const std::string&                   well_name,
                                const std::vector<data::Rates::opt>& phs,
                                const data::Well&                    rst_well,
                                SingleWellState&                     ws) const;

    //! \brief Loads per-well data from restart structures.
    void loadRestartWellData(const std::string&                   well_name,
                             const bool                           handle_ms_well,
                             const std::vector<data::Rates::opt>& phs,
                             const data::Well&                    rst_well,
                             const std::vector<PerforationData>&  old_perf_data,
                             SingleWellState&                     ws) const;

    //! \brief Loads per-group data from restart structures.
    void loadRestartGroupData(const std::string&     group,
                              const data::GroupData& value,
                              GroupState& grpState) const;

    const BlackoilWellModelGeneric& wellModel_; //!< Reference to well model
};


} // namespace Opm

#endif
