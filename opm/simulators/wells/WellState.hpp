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
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm
{
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
                  const SummaryState& summary_state)
        {
            // clear old name mapping
            wellMap_.clear();

            well_perf_data_ = well_perf_data;
            parallel_well_info_ = parallel_well_info;

            {
                // const int nw = wells->number_of_wells;
                const int nw = wells_ecl.size();
                const int np = this->phase_usage_.num_phases;
                // const int np = wells->number_of_phases;
                status_.assign(nw, Well::Status::OPEN);
                bhp_.resize(nw, 0.0);
                thp_.resize(nw, 0.0);
                temperature_.resize(nw, 273.15 + 15.56); // standard condition temperature
                wellrates_.resize(nw * np, 0.0);
                int connpos = 0;
                for (int w = 0; w < nw; ++w) {
                    const Well& well = wells_ecl[w];
              
                    // Initialize bhp(), thp(), wellRates(), temperature().
                    initSingleWell(cellPressures, w, well, *parallel_well_info[w], summary_state);

                    // Setup wellname -> well index mapping.
                    const int num_perf_this_well = well_perf_data[w].size();
                    std::string name = well.name();
                    assert( name.size() > 0 );
                    mapentry_t& wellMapEntry = wellMap_[name];
                    wellMapEntry[ 0 ] = w;
                    wellMapEntry[ 1 ] = connpos;
                    wellMapEntry[ 2 ] = num_perf_this_well;
                    connpos += num_perf_this_well;
                }

                // The perforation rates and perforation pressures are
                // not expected to be consistent with bhp_ and wellrates_
                // after init().
                perfrates_.resize(connpos, 0.0);
                perfpress_.resize(connpos, -1e100);
            }
        }

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
                                         const std::vector<PerforationData>& well_perf_data)
        {
            if (this->well_perf_data_[well_index].size() != well_perf_data.size()) {
                throw std::invalid_argument {
                    "Size mismatch for perforation data in well "
                    + std::to_string(well_index)
                };
            }

            auto connID = std::size_t{0};
            auto dst = this->well_perf_data_[well_index].begin();
            for (const auto& src : well_perf_data) {
                if (dst->cell_index != src.cell_index) {
                    throw std::invalid_argument {
                        "Cell index mismatch in connection "
                        + std::to_string(connID)
                        + " of well "
                        + std::to_string(well_index)
                    };
                }

                if (dst->satnum_id != src.satnum_id) {
                    throw std::invalid_argument {
                        "Saturation function table mismatch in connection "
                        + std::to_string(connID)
                        + " of well "
                        + std::to_string(well_index)
                    };
                }

                dst->connection_transmissibility_factor =
                    src.connection_transmissibility_factor;

                ++dst;
                ++connID;
            }
        }

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

        const ParallelWellInfo& parallelWellInfo(std::size_t well_index) const
        {
            return *parallel_well_info_[well_index];
        }

        bool wellIsOwned(std::size_t well_index, [[maybe_unused]] const std::string& wellName) const
        {
            const auto& well_info = parallelWellInfo(well_index);
            assert(well_info.name() == wellName);

            return well_info.isOwner();
        }

        bool wellIsOwned(const std::string& wellName) const
        {
            const auto& it = wellMap().find( wellName );
            if (it == wellMap().end()) {
                OPM_THROW(std::logic_error, "Could not find well " << wellName << " in well map");
            }
            const int well_index = it->second[0];
            return wellIsOwned(well_index, wellName);
        }

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

        virtual void shutWell(int well_index) {
            this->status_[well_index] = Well::Status::SHUT;
            this->thp_[well_index] = 0;
            this->bhp_[well_index] = 0;
            const int np = numPhases();
            for (int p = 0; p < np; ++p)
                this->wellrates_[np * well_index + p] = 0;
        }

        virtual void stopWell(int well_index) {
            this->status_[well_index] = Well::Status::STOP;
            this->thp_[well_index] = 0;
        }

        void updateStatus(int well_index, Well::Status status) {
            switch (status) {
            case Well::Status::OPEN:
                this->openWell(well_index);
                break;
            case Well::Status::SHUT:
                this->shutWell(well_index);
                break;
            case Well::Status::STOP:
                this->stopWell(well_index);
                break;
            default:
                throw std::logic_error("Invalid well status");
            }
        }

        virtual data::Wells
        report(const int* globalCellIdxMap,
               const std::function<bool(const int)>& wasDynamicallyClosed) const
        {
            using rt = data::Rates::opt;

            const auto& pu = this->phaseUsage();
            data::Wells dw;
            for( const auto& itr : this->wellMap_ ) {
                const auto well_index = itr.second[ 0 ];
                if ((this->status_[well_index] == Well::Status::SHUT) &&
                    ! wasDynamicallyClosed(well_index))
                {
                    continue;
                }

                const auto& pwinfo = *parallel_well_info_[well_index];
                using WellT = std::remove_reference_t<decltype(dw[ itr.first ])>;
                WellT dummyWell; // dummy if we are not owner
                auto& well = pwinfo.isOwner() ? dw[ itr.first ] : dummyWell;
                well.bhp = this->bhp().at( well_index );
                well.thp = this->thp().at( well_index );
                well.temperature = this->temperature().at( well_index );

                const auto wellrate_index = well_index * pu.num_phases;
                const auto& wv = this->wellRates();
                if( pu.phase_used[BlackoilPhases::Aqua] ) {
                    well.rates.set( rt::wat, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Aqua] ] );
                }

                if( pu.phase_used[BlackoilPhases::Liquid] ) {
                    well.rates.set( rt::oil, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Liquid] ] );
                }

                if( pu.phase_used[BlackoilPhases::Vapour] ) {
                    well.rates.set( rt::gas, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Vapour] ] );
                }

                if (pwinfo.communication().size()==1)
                {
                    reportConnections(well, pu, itr, globalCellIdxMap);
                }
                else
                {
                    assert(pwinfo.communication().rank() != 0 || &dummyWell != &well);
                    // report the local connections
                    reportConnections(dummyWell, pu, itr, globalCellIdxMap);
                    // gather them to well on root.
                    gatherVectorsOnRoot(dummyWell.connections, well.connections,
                                        pwinfo.communication());
                }
            }

            return dw;

        }

        virtual void reportConnections(data::Well& well, [[maybe_unused]] const PhaseUsage& pu,
                                       const WellMapType::value_type& itr,
                                       const int* globalCellIdxMap) const
        {
            const auto well_index = itr.second[ 0 ];
            const auto& pd = this->well_perf_data_[well_index];
            const int num_perf_well = pd.size();
            well.connections.resize(num_perf_well);

            const auto * perf_rates = &this->perfRates()[itr.second[1]];
            const auto * perf_pressure = &this->perfPress()[itr.second[1]];
            for( int i = 0; i < num_perf_well; ++i ) {
                const auto active_index = this->well_perf_data_[well_index][i].cell_index;
                auto& connection = well.connections[ i ];
                connection.index = globalCellIdxMap[active_index];
                connection.pressure = perf_pressure[i];
                connection.reservoir_rate = perf_rates[i];
                connection.trans_factor = pd[i].connection_transmissibility_factor;
            }
            assert(num_perf_well == int(well.connections.size()));
        }
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

        using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
        using Communication = Dune::Communication<MPIComm>;
#else
        using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif
        void gatherVectorsOnRoot(const std::vector< data::Connection >& from_connections,
                                 std::vector< data::Connection >& to_connections,
                                 const Communication& comm) const
        {
            int size = from_connections.size();
            std::vector<int> sizes;
            std::vector<int> displ;
            if (comm.rank()==0){
                sizes.resize(comm.size());
            }
            comm.gather(&size, sizes.data(), 1, 0);

            if (comm.rank()==0){
                displ.resize(comm.size()+1, 0);
                std::partial_sum(sizes.begin(), sizes.end(), displ.begin()+1);
                to_connections.resize(displ.back());
            }
            comm.gatherv(from_connections.data(), size, to_connections.data(),
                         sizes.data(), displ.data(), 0);
        }
        void initSingleWell(const std::vector<double>& cellPressures,
                            const int w,
                            const Well& well,
                            const ParallelWellInfo& well_info,
                            const SummaryState& summary_state)
        {
            assert(well.isInjector() || well.isProducer());

            // Set default zero initial well rates.
            // May be overwritten below.
            const auto& pu = this->phase_usage_;
            const int np = pu.num_phases;
            for (int p = 0; p < np; ++p) {
                wellrates_[np*w + p] = 0.0;
            }

            if ( well.isInjector() ) { 
                temperature_[w] = well.injectionControls(summary_state).temperature;
            }

            const int num_perf_this_well = well_info.communication().sum(well_perf_data_[w].size());
            if ( num_perf_this_well == 0 ) {
                // No perforations of the well. Initialize to zero.
                bhp_[w] = 0.;
                thp_[w] = 0.;
                return;
            }

            const auto inj_controls = well.isInjector() ? well.injectionControls(summary_state) : Well::InjectionControls(0);
            const auto prod_controls = well.isProducer() ? well.productionControls(summary_state) : Well::ProductionControls(0);

            const bool is_bhp = well.isInjector() ? (inj_controls.cmode == Well::InjectorCMode::BHP)
                : (prod_controls.cmode == Well::ProducerCMode::BHP);
            const double bhp_limit = well.isInjector() ? inj_controls.bhp_limit : prod_controls.bhp_limit;
            const bool is_grup = well.isInjector() ? (inj_controls.cmode == Well::InjectorCMode::GRUP)
                : (prod_controls.cmode == Well::ProducerCMode::GRUP);

            const double inj_surf_rate = well.isInjector() ? inj_controls.surface_rate : 0.0; // To avoid a "maybe-uninitialized" warning.

            const double local_pressure = well_perf_data_[w].empty() ?
                0 : cellPressures[well_perf_data_[w][0].cell_index];
            const double global_pressure = well_info.broadcastFirstPerforationValue(local_pressure);

            if (well.getStatus() == Well::Status::OPEN) {
                this->openWell(w);
            }

            if (well.getStatus() == Well::Status::STOP) {
                // Stopped well:
                // 1. Rates: zero well rates.
                // 2. Bhp: assign bhp equal to bhp control, if
                //    applicable, otherwise assign equal to
                //    first perforation cell pressure.
                if (is_bhp) {
                    bhp_[w] = bhp_limit;
                } else {
                    bhp_[w] = global_pressure;
                }
            } else if (is_grup) {
                // Well under group control.
                // 1. Rates: zero well rates.
                // 2. Bhp: initialize bhp to be a
                //    little above or below (depending on if
                //    the well is an injector or producer)
                //    pressure in first perforation cell.
                const double safety_factor = well.isInjector() ? 1.01 : 0.99;
                bhp_[w] = safety_factor * global_pressure;
            } else {
                // Open well, under own control:
                // 1. Rates: initialize well rates to match
                //    controls if type is ORAT/GRAT/WRAT
                //    (producer) or RATE (injector).
                //    Otherwise, we cannot set the correct
                //    value here and initialize to zero rate.
                if (well.isInjector()) {
                    if (inj_controls.cmode == Well::InjectorCMode::RATE) {
                        switch (inj_controls.injector_type) {
                        case InjectorType::WATER:
                            assert(pu.phase_used[BlackoilPhases::Aqua]);
                            wellrates_[np*w + pu.phase_pos[BlackoilPhases::Aqua]] = inj_surf_rate;
                            break;
                        case InjectorType::GAS:
                            assert(pu.phase_used[BlackoilPhases::Vapour]);
                            wellrates_[np*w + pu.phase_pos[BlackoilPhases::Vapour]] = inj_surf_rate;
                            break;
                        case InjectorType::OIL:
                            assert(pu.phase_used[BlackoilPhases::Liquid]);
                            wellrates_[np*w + pu.phase_pos[BlackoilPhases::Liquid]] = inj_surf_rate;
                            break;
                        case InjectorType::MULTI:
                            // Not currently handled, keep zero init.
                            break;
                        }
                    } else {
                        // Keep zero init.
                    }
                } else {
                    assert(well.isProducer());
                    // Note negative rates for producing wells.
                    switch (prod_controls.cmode) {
                    case Well::ProducerCMode::ORAT:
                        assert(pu.phase_used[BlackoilPhases::Liquid]);
                        wellrates_[np*w + pu.phase_pos[BlackoilPhases::Liquid]] = -prod_controls.oil_rate;
                        break;
                    case Well::ProducerCMode::WRAT:
                        assert(pu.phase_used[BlackoilPhases::Aqua]);
                        wellrates_[np*w + pu.phase_pos[BlackoilPhases::Aqua]] = -prod_controls.water_rate;
                        break;
                    case Well::ProducerCMode::GRAT:
                        assert(pu.phase_used[BlackoilPhases::Vapour]);
                        wellrates_[np*w + pu.phase_pos[BlackoilPhases::Vapour]] = -prod_controls.gas_rate;
                        break;
                    default:
                        // Keep zero init.
                        break;
                    }
                }
                // 2. Bhp: initialize bhp to be target pressure if
                //    bhp-controlled well, otherwise set to a
                //    little above or below (depending on if
                //    the well is an injector or producer)
                //    pressure in first perforation cell.
                if (is_bhp) {
                    bhp_[w] = bhp_limit;
                } else {
                    const double safety_factor = well.isInjector() ? 1.01 : 0.99;
                    bhp_[w] = safety_factor * global_pressure;
                }
            }

            // 3. Thp: assign thp equal to thp target/limit, if such a limit exists,
            //    otherwise keep it zero.
            const bool has_thp = well.isInjector() ? inj_controls.hasControl(Well::InjectorCMode::THP)
                : prod_controls.hasControl(Well::ProducerCMode::THP);
            const double thp_limit = well.isInjector() ? inj_controls.thp_limit : prod_controls.thp_limit;
            if (has_thp) {
                thp_[w] = thp_limit;
            }

        }

    protected:
        std::vector<std::vector<PerforationData>> well_perf_data_;
        std::vector<ParallelWellInfo*> parallel_well_info_;
    };
    
} // namespace Opm

#endif // OPM_WELLSTATE_HEADER_INCLUDED
