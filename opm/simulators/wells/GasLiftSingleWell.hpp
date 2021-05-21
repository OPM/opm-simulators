/*
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
#define OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED

// NOTE: StandardWell.hpp includes ourself (GasLiftSingleWell.hpp), so we need
//   to forward declare StandardWell for it to be defined in this file.
namespace Opm {
    template<typename TypeTag> class StandardWell;
}
#include <opm/simulators/wells/StandardWell.hpp>

#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>

#include <optional>
#include <vector>
#include <utility>
#include <fmt/format.h>

namespace Opm
{
    template<class TypeTag>
    class GasLiftSingleWell : public GasLiftSingleWellGeneric
    {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using WellState = WellStateFullyImplicitBlackoil;
        using StdWell = StandardWell<TypeTag>;

    public:
        GasLiftSingleWell(
            const StdWell &std_well,
            const Simulator &ebos_simulator,
            const SummaryState &summary_state,
            DeferredLogger &deferred_logger,
            WellState &well_state
        );
        const WellInterface<TypeTag> &getStdWell() const { return std_well_; }

    private:
        std::optional<double> computeBhpAtThpLimit_(double alq) const override;
        void computeWellRates_(
            double bhp, std::vector<double> &potentials, bool debug_output=true) const override;

        void setAlqMaxRate_(const GasLiftOpt::Well& well);

        const Simulator &ebos_simulator_;
        const StdWell &std_well_;
    };

} // namespace Opm

#include "GasLiftSingleWell_impl.hpp"

#endif // OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
