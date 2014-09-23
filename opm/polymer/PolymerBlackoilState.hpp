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

#ifndef OPM_POLYMERBLACKOILSTATE_HEADER_INCLUDED
#define OPM_POLYMERBLACKOILSTATE_HEADER_INCLUDED


#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/grid.h>
#include <vector>

namespace Opm
{

    /// Simulator state for a compressible two-phase simulator with polymer.
    /// We use the Blackoil state parameters.
    class PolymerBlackoilState
    {
    public:
        void init(const UnstructuredGrid& g, int num_phases)
        {
            this->init(g.number_of_cells, g.number_of_faces, num_phases);
        }

        void init(int number_of_cells, int number_of_faces, int num_phases)
        {
            state_blackoil_.init(number_of_cells, number_of_faces, num_phases);
            concentration_.resize(number_of_cells, 0.0);
            cmax_.resize(number_of_cells, 0.0);
        }
        int numPhases() const
        {
            return state_blackoil_.numPhases();
        }

        enum ExtremalSat { MinSat = BlackoilState::MinSat, MaxSat = BlackoilState::MaxSat };

        void setFirstSat(const std::vector<int>& cells,
                         const Opm::BlackoilPropertiesInterface& props,
                         ExtremalSat es)
        {
            // A better solution for embedding BlackoilState::ExtremalSat could perhaps
            // be found, to avoid the cast.
            state_blackoil_.setFirstSat(cells, props, static_cast<BlackoilState::ExtremalSat>(es));
        }

        std::vector<double>& pressure    ()     { return state_blackoil_.pressure(); }
        std::vector<double>& surfacevol  ()     { return state_blackoil_.surfacevol(); }
        std::vector<double>& facepressure()     { return state_blackoil_.facepressure(); }
        std::vector<double>& faceflux    ()     { return state_blackoil_.faceflux(); }
        std::vector<double>& saturation  ()     { return state_blackoil_.saturation(); }
        std::vector<double>& concentration()    { return concentration_; }
        std::vector<double>& maxconcentration() { return cmax_; }

        const std::vector<double>& pressure    () const     { return state_blackoil_.pressure(); }
        const std::vector<double>& surfacevol  () const     { return state_blackoil_.surfacevol(); }
        const std::vector<double>& facepressure() const     { return state_blackoil_.facepressure(); }
        const std::vector<double>& faceflux    () const     { return state_blackoil_.faceflux(); }
        const std::vector<double>& saturation  () const     { return state_blackoil_.saturation(); }
        const std::vector<double>& concentration() const    { return concentration_; }
        const std::vector<double>& maxconcentration() const { return cmax_; }

        BlackoilState& blackoilState() { return state_blackoil_; }
        const BlackoilState& blackoilState() const { return state_blackoil_; }

    private:
        BlackoilState state_blackoil_;
        std::vector<double> concentration_;
        std::vector<double> cmax_;
    };

} // namespace Opm




#endif // OPM_POLYMERBLACKOILSTATE_HEADER_INCLUDED
