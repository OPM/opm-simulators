/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELLSTATE_HEADER_INCLUDED
#define OPM_WELLSTATE_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/wells/PerforationData.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace Opm
{
    class ParallelWellInfo;
    class SummaryState;

    /// The state of a set of wells.
    class WellState
    {
    public:
        typedef std::array< int, 3 >  mapentry_t;
        typedef std::map< std::string, mapentry_t > WellMapType;



        explicit WellState(const PhaseUsage& pu) :
            phase_usage_(pu)
        {}


        /// Allocate and initialize if wells is non-null.
        /// Also tries to give useful initial values to the bhp() and
        /// wellRates() fields, depending on controls.  The
        /// perfRates() field is filled with zero, and perfPress()
        /// with -1e100.
        void init(const std::vector<double>& cellPressures,
                  const std::vector<Well>& wells_ecl,
                  const std::vector<ParallelWellInfo*>& parallel_well_info,
                  const std::vector<std::vector<PerforationData>>& well_perf_data,
                  const SummaryState& summary_state);

        /// Special purpose method to support dynamically rescaling a well's
        /// CTFs through WELPI.
        ///
        /// \param[in] well_index Process-local linear index of single well.
        ///    Must be in the range 0..numWells()-1.
        ///
        /// \param[in] well_perf_data New perforation data.  Only
        ///    PerforationData::connection_transmissibility_factor actually
        ///    used (overwrites existing internal values).
        void resetConnectionTransFactors(const int well_index,
                                         const std::vector<PerforationData>& well_perf_data);

        /// One bhp pressure per well.
        std::vector<double>& bhp() { return bhp_; }
        const std::vector<double>& bhp() const { return bhp_; }

        /// One thp pressure per well.
        std::vector<double>& thp() { return thp_; }
        const std::vector<double>& thp() const { return thp_; }

        /// One temperature per well.
        std::vector<double>& temperature() { return temperature_; }
        const std::vector<double>& temperature() const { return temperature_; }

        /// One rate per well and phase.
        std::vector<double>& wellRates() { return wellrates_; }
        const std::vector<double>& wellRates() const { return wellrates_; }

        /// One rate per well connection.
        std::vector<double>& perfRates() { return perfrates_; }
        const std::vector<double>& perfRates() const { return perfrates_; }

        /// One pressure per well connection.
        std::vector<double>& perfPress() { return perfpress_; }
        const std::vector<double>& perfPress() const { return perfpress_; }

        const WellMapType& wellMap() const { return wellMap_; }
        WellMapType& wellMap() { return wellMap_; }

        const ParallelWellInfo& parallelWellInfo(std::size_t well_index) const;

        bool wellIsOwned(std::size_t well_index,
                         const std::string& wellName) const;

        bool wellIsOwned(const std::string& wellName) const;

        /// The number of wells present.
        int numWells() const
        {
            return bhp().size();
        }

        /// The number of phases present.
        int numPhases() const
        {
            return this->phase_usage_.num_phases;
        }

        const PhaseUsage& phaseUsage() const {
            return this->phase_usage_;
        }



        void openWell(int well_index) {
            this->status_[well_index] = Well::Status::OPEN;
        }

        virtual void shutWell(int well_index);

        virtual void stopWell(int well_index);

        void updateStatus(int well_index, Well::Status status);

        virtual data::Wells
        report(const int* globalCellIdxMap,
               const std::function<bool(const int)>& wasDynamicallyClosed) const;

        virtual void reportConnections(data::Well& well, const PhaseUsage&,
                                       const WellMapType::value_type& itr,
                                       const int* globalCellIdxMap) const;
        virtual ~WellState() = default;
        WellState() = default;
        WellState(const WellState& rhs)  = default;
        WellState& operator=(const WellState& rhs) = default;

    private:
        PhaseUsage phase_usage_;
        std::vector<double> bhp_;
        std::vector<double> thp_;
        std::vector<double> temperature_;
        std::vector<double> wellrates_;
        std::vector<double> perfrates_;
        std::vector<double> perfpress_;
    protected:
        std::vector<Well::Status> status_;
    private:

        WellMapType wellMap_;

        template<class Communication>
        void gatherVectorsOnRoot(const std::vector< data::Connection >& from_connections,
                                 std::vector< data::Connection >& to_connections,
                                 const Communication& comm) const;

        void initSingleWell(const std::vector<double>& cellPressures,
                            const int w,
                            const Well& well,
                            const ParallelWellInfo& well_info,
                            const SummaryState& summary_state);

    protected:
        std::vector<std::vector<PerforationData>> well_perf_data_;
        std::vector<ParallelWellInfo*> parallel_well_info_;
    };

} // namespace Opm

#endif // OPM_WELLSTATE_HEADER_INCLUDED
