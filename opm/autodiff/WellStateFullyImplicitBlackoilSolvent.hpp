/*
  Copyright 2015 IRIS AS

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

#ifndef OPM_WELLSTATEFULLYIMPLICITBLACKOILSOLVENT_HEADER_INCLUDED
#define OPM_WELLSTATEFULLYIMPLICITBLACKOILSOLVENT_HEADER_INCLUDED

#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

namespace Opm
{

    class WellStateFullyImplicitBlackoilSolvent : public WellStateFullyImplicitBlackoil
    {
        typedef WellStateFullyImplicitBlackoil  BaseType;
    public:

        /// One solvent fraction per well connection
        std::vector<double>& solventFraction() { return solvent_fraction_; }
        const std::vector<double>& solventFraction() const { return solvent_fraction_; }

        data::Wells report(const PhaseUsage &pu) const override {
            data::Wells res = WellStateFullyImplicitBlackoil::report(pu);

            const int nw = WellState::numWells();

            // If there are now wells numPhases throws a floating point
            // exception.
            if (nw == 0) {
                return res;
            }

            const int np = BaseType::numPhases();

            assert( np == 3 ); // the solvent model assumes 3 phases in the base model

            // completions aren't supported yet
            for( auto w = 0; w < nw; ++w ) {
                using rt = data::Rates::opt;
                double solvent_well_rate = 0.0;
                for (int perf = wells_->well_connpos[w]; perf < wells_->well_connpos[w+1]; ++perf ) {
                    auto solvent_rate_this = BaseType::perfPhaseRates()[np*perf + pu.phase_pos[BlackoilPhases::Vapour]] * solventFraction()[perf];
                    solvent_well_rate += solvent_rate_this;
                }

                res.at( wells_->name[ w ]).rates.set( rt::solvent, solvent_well_rate );

            }

            return res;
        }
    private:
        std::vector<double> solvent_fraction_;
    };

} // namespace Opm

#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOILSOLVENT_HEADER_INCLUDED
