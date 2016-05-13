/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
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

#ifndef OPM_BLACKOILSTATE_HEADER_INCLUDED
#define OPM_BLACKOILSTATE_HEADER_INCLUDED

#include <opm/common/data/SimulationDataContainer.hpp>

#include <opm/core/grid.h>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <vector>

namespace Opm
{

    enum HydroCarbonState {
        GasOnly = 0,
        GasAndOil = 1,
        OilOnly = 2
    };

    /// Simulator state for a blackoil simulator.
    class BlackoilState : public SimulationDataContainer
    {
    public:
        static const std::string GASOILRATIO;
        static const std::string RV;
        static const std::string SURFACEVOL;

        /// Main constructor setting the sizes for the contained data
        /// types.
        /// \param num_cells   number of elements in cell data vectors
        /// \param num_faces   number of elements in face data vectors
        /// \param num_phases  number of phases, the number of components
        ///                    in any data vector must equal 1 or this
        ///                    number (this behaviour and argument is deprecated).
        BlackoilState(size_t num_cells, size_t num_faces, size_t num_phases);

        /// Copy constructor.
        /// Must be defined explicitly because class contains non-value objects
        /// (the reference pointers rv_ref_ etc.) that should not simply
        /// be copied.
        BlackoilState(const BlackoilState& other);

        /// Copy assignment operator.
        /// Must be defined explicitly because class contains non-value objects
        /// (the reference pointers rv_ref_ etc.) that should not simply
        /// be copied.
        BlackoilState& operator=(const BlackoilState& other);

        std::vector<double>& surfacevol  () { return *surfacevol_ref_;  }
        std::vector<double>& gasoilratio () { return *gasoilratio_ref_; }
        std::vector<double>& rv ()          { return *rv_ref_;          }
        std::vector<HydroCarbonState>& hydroCarbonState() { return hydrocarbonstate_;  }

        const std::vector<double>& surfacevol  () const { return *surfacevol_ref_;  }
        const std::vector<double>& gasoilratio () const { return *gasoilratio_ref_; }
        const std::vector<double>& rv ()          const { return *rv_ref_;          }
        const std::vector<HydroCarbonState>& hydroCarbonState() const { return hydrocarbonstate_;  }

    private:
        void setBlackoilStateReferencePointers();
        std::vector<double>* surfacevol_ref_;
        std::vector<double>* gasoilratio_ref_;
        std::vector<double>* rv_ref_;

        // A vector storing the hydro carbon state.
        std::vector<HydroCarbonState> hydrocarbonstate_;


    };
} // namespace Opm


#endif // OPM_BLACKOILSTATE_HEADER_INCLUDED
