/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

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


#ifndef OPM_STANDARDWELLSSOLVENT_HEADER_INCLUDED
#define OPM_STANDARDWELLSSOLVENT_HEADER_INCLUDED

#include <opm/autodiff/StandardWells.hpp>
#include <opm/autodiff/SolventPropsAdFromDeck.hpp>

namespace Opm {


        /// Class for handling the standard well model for solvent model
        class StandardWellsSolvent : public StandardWells
        {
        public:

            using Base = StandardWells;
            using Base::computeWellConnectionDensitesPressures;

            // ---------  Public methods  ---------
            StandardWellsSolvent(const Wells* wells_arg);

            // added the Solvent related
            void initSolvent(const SolventPropsAdFromDeck* solvent_props,
                             const int solvent_pos,
                             const bool has_solvent);

            template <class SolutionState, class WellState>
            void computePropertiesForWellConnectionPressures(const SolutionState& state,
                                                             const WellState& xw,
                                                             std::vector<double>& b_perf,
                                                             std::vector<double>& rsmax_perf,
                                                             std::vector<double>& rvmax_perf,
                                                             std::vector<double>& surf_dens_perf);

            // TODO: fluid and active may be can put in the member list
            template <class ReservoirResidualQuant, class SolutionState>
            void extractWellPerfProperties(const SolutionState& state,
                                           const std::vector<ReservoirResidualQuant>& rq,
                                           std::vector<ADB>& mob_perfcells,
                                           std::vector<ADB>& b_perfcells) const;

            template <class SolutionState, class WellState>
            void computeWellConnectionPressures(const SolutionState& state,
                                                const WellState& xw,
                                                const Vector& depth,
                                                const double gravity);

        protected:
            const SolventPropsAdFromDeck* solvent_props_;
            int solvent_pos_;
            bool has_solvent_;

            using Base::phase_condition_;

        };


} // namespace Opm

#include "StandardWellsSolvent_impl.hpp"

#endif
