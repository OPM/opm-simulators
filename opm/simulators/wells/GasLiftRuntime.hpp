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

#ifndef OPM_GASLIFT_RUNTIME_HEADER_INCLUDED
#define OPM_GASLIFT_RUNTIME_HEADER_INCLUDED

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
// NOTE: StandardWell.hpp includes ourself (GasLiftRuntime.hpp), so we need
//   to forward declare StandardWell for it to be defined in this file.
namespace Opm {
    template<typename TypeTag> class StandardWell;
}
#include <opm/simulators/wells/StandardWell.hpp>

#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <cassert>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <fmt/format.h>

namespace Opm
{
    template<class TypeTag>
    class GasLiftRuntime {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using WellState = WellStateFullyImplicitBlackoil;
        using StdWell = Opm::StandardWell<TypeTag>;
        // TODO: same definition with WellInterface, and
        //  WellStateFullyImplicitBlackoil eventually they should go
        //  to a common header file.
        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;
        struct OptimizeState;
    public:
        GasLiftRuntime(
            const StdWell &std_well,
            const Simulator &ebos_simulator,
            const SummaryState &summary_state,
            DeferredLogger &deferred_logger,
            const WellState &well_state,
            const Well::ProductionControls &controls
        );
        void runOptimize();
    private:
        void computeInitialWellRates_();
        void computeWellRates_(double bhp, std::vector<double> &potentials);
        void debugShowIterationInfo_(OptimizeState &state, double alq);
        void debugShowStartIteration_(double alq, bool increase);
        void displayDebugMessage_(const std::string &msg);
        void displayWarning_();
        void displayWarning_(std::string warning);
        double getGasRateWithLimit_(std::vector<double> &potentials);
        double getOilRateWithLimit_(std::vector<double> &potentials);
        void logSuccess_();
        bool runOptimizeLoop_(bool increase);
        void setAlqMaxRate_(const GasLiftOpt::Well &well);
        void setAlqMinRate_(const GasLiftOpt::Well &well);
        bool tryIncreaseLiftGas_();
        bool tryDecreaseLiftGas_();
        void updateWellStateAlqFixedValue_(const GasLiftOpt::Well &well);
        bool useFixedAlq_(const GasLiftOpt::Well &well);
        void warnMaxIterationsExceeded_();

        const Well::ProductionControls &controls_;
        DeferredLogger &deferred_logger_;
        const Simulator &ebos_simulator_;
        std::vector<double> potentials_;
        const StdWell &std_well_;
        const SummaryState &summary_state_;
        const WellState &well_state_;
        std::string well_name_;
        bool debug;  // extra debug output

        double alpha_w_;
        double alpha_g_;
        double eco_grad_;
        int gas_pos_;
        bool has_run_init_ = false;
        double increment_;
        double max_alq_;
        int max_iterations_;
        double min_alq_;
        double new_alq_;
        int oil_pos_;
        bool optimize_;
        double orig_alq_;
        int water_pos_;

        struct OptimizeState {
            OptimizeState( GasLiftRuntime &parent_, bool increase_ ) :
                parent(parent_),
                increase(increase_),
                it(0),
                stop_iteration(false),
                bhp(-1)
            {}

            GasLiftRuntime &parent;
            bool increase;
            int it;
            bool stop_iteration;
            double bhp;

            double addOrSubtractAlqIncrement(double alq);
            double calcGradient(double oil_rate, double new_oil_rate,
                double gas_rate, double new_gas_rate);
            bool checkAlqOutsideLimits(double alq, double oil_rate);
            bool checkEcoGradient(double gradient);
            bool checkOilRateExceedsTarget(double oil_rate);
            bool checkRate(double rate, double limit, const std::string rate_str);
            bool checkWellRatesViolated(std::vector<double> &potentials);
            bool computeBhpAtThpLimit(double alq);
            double getBhpWithLimit();
            void warn_(std::string msg) {parent.displayWarning_(msg);}
        };

    };

} // namespace Opm

#include "GasLiftRuntime_impl.hpp"

#endif // OPM_GASLIFT_RUNTIME_HEADER_INCLUDED
