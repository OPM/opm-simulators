/*
  Copyright 2021 Equinor ASA.

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

#ifndef OPM_PERFDATA_HEADER_INCLUDED
#define OPM_PERFDATA_HEADER_INCLUDED

#include <opm/simulators/wells/ConnFiltrateData.hpp>
#include <opm/simulators/wells/ConnFracStatistics.hpp>

#include <cstddef>
#include <vector>
#include <array>

namespace Opm {

template<class Scalar>
class PerfData
{
private:
    bool injector{false};

public:
    PerfData() = default;
    PerfData(std::size_t num_perf,
             Scalar pressure_first_connection_,
             bool injector_,
             std::size_t num_phases);

    static PerfData serializationTestObject();

    std::size_t size() const;
    bool empty() const;
    bool try_assign(const PerfData& other);

    /// \brief Make containers valid for injectors
    ///
    /// Needed if a producer is switched to an injector,
    void prepareInjectorContainers();

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(injector);
        serializer(pressure_first_connection);
        serializer(pressure);
        serializer(rates);
        serializer(phase_rates);
        serializer(phase_mixing_rates);
        serializer(solvent_rates);
        serializer(polymer_rates);
        serializer(brine_rates);
        serializer(prod_index);
        serializer(microbial_rates);
        serializer(oxygen_rates);
        serializer(urea_rates);
        serializer(cell_index);
        serializer(connection_transmissibility_factor);
        serializer(connection_d_factor);
        serializer(connection_compaction_tmult);
        serializer(satnum_id);
        serializer(ecl_index);
        serializer(water_throughput);
        serializer(skin_pressure);
        serializer(water_velocity);
        serializer(filtrate_data);
        serializer(connFracStatistics);
        serializer(gas_mass_rates);
    }

    bool operator==(const PerfData&) const;

    // Note to maintainers: If you make changes to this list of data
    // members, then please update the constructor, operator==(),
    // serializationTestObject(), and serializeOp() accordingly.  Moreover,
    // if you're adding a new member representing a dynamically calculated
    // result, e.g., a flow rate, then please update try_assign() as well.

    Scalar pressure_first_connection{};
    std::vector<Scalar> pressure{};
    std::vector<Scalar> rates{};
    std::vector<Scalar> phase_rates{};
    std::vector<std::array<Scalar,4>> phase_mixing_rates{};
    std::vector<Scalar> solvent_rates{};
    std::vector<Scalar> polymer_rates{};
    std::vector<Scalar> brine_rates{};
    std::vector<Scalar> prod_index{};
    std::vector<Scalar> microbial_rates{};
    std::vector<Scalar> oxygen_rates{};
    std::vector<Scalar> urea_rates{};
    std::vector<std::size_t> cell_index{};
    std::vector<Scalar> connection_transmissibility_factor{};
    std::vector<Scalar> connection_d_factor{};
    std::vector<Scalar> connection_compaction_tmult{};
    std::vector<int> satnum_id{};
    std::vector<std::size_t> ecl_index{};
    std::vector<Scalar> gas_mass_rates{};

    // The water_throughput, skin_pressure and water_velocity variables are
    // only used for injectors to check the injectivity.
    std::vector<Scalar> water_throughput{};
    std::vector<Scalar> skin_pressure{};
    std::vector<Scalar> water_velocity{};

    ConnFiltrateData<Scalar> filtrate_data{};
    std::vector<ConnFracStatistics<Scalar>> connFracStatistics{};
};

} // namespace Opm

#endif // OPM_PERFDATA_HEADER_INCLUDED
