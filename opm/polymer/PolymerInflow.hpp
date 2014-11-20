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

#ifndef OPM_POLYMERINFLOW_HEADER_INCLUDED
#define OPM_POLYMERINFLOW_HEADER_INCLUDED

#include <opm/core/utility/SparseVector.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/WellPolymerProperties.hpp>
#include <vector>
#include <string>
#include <unordered_map>

struct Wells;

namespace Opm
{
    /// @brief Interface for classes encapsulating polymer inflow information.
    class PolymerInflowInterface
    {
    public:
        /// Virtual destructor for subclassing.
        virtual ~PolymerInflowInterface() {}

        /// Get inflow concentrations for all cells.
        /// \param[in]  step_start     Start of timestep.
        /// \param[in]  step_end       End of timestep.
        /// \param[out] poly_inflow_c  Injection concentrations to use for timestep, per cell.
        ///                            Must be properly sized before calling.
        virtual void getInflowValues(const double step_start,
                                     const double step_end,
                                     std::vector<double>& poly_inflow_c) const = 0;
    };



    /// @brief Basic polymer injection behaviour class.
    /// This class gives all injectors the same polymer concentration,
    /// during a single time interval. Amount and interval can be specified.
    class PolymerInflowBasic : public PolymerInflowInterface
    {
    public:
        /// Constructor.
        /// \param[in]  starttime  Start time of injection in seconds.
        /// \param[in]  endtime    End time of injection in seconds.
        /// \param[in]  amount     Amount to be injected per second.
        PolymerInflowBasic(const double starttime,
                           const double endtime,
                           const double amount);

        /// Get inflow concentrations for all cells.
        /// \param[in]  step_start     Start of timestep.
        /// \param[in]  step_end       End of timestep.
        /// \param[out] poly_inflow_c  Injection concentrations to use for timestep, per cell.
        ///                            Must be properly sized before calling.
        virtual void getInflowValues(const double step_start,
                                     const double step_end,
                                     std::vector<double>& poly_inflow_c) const;
    private:
        double stime_;
        double etime_;
        double amount_;
    };


    /// @brief Polymer injection behaviour class using deck WPOLYMER.
    /// This class reads the accumulated WPOLYMER lines from the deck,
    /// and applies the last row given for each well.
    class PolymerInflowFromDeck : public PolymerInflowInterface
    {
    public:
        /// Constructor.
        /// \param[in]  deck        Input deck expected to contain WPOLYMER.
        /// \param[in]  wells       Wells structure.
        /// \param[in]  num_cells   Number of cells in grid.
        PolymerInflowFromDeck(Opm::DeckConstPtr deck,
                              const Wells& wells,
                              const int num_cells);

        /// Constructor.
        /// \param[in]  deck        Input deck expected to contain WPOLYMER.
        /// \param[in]  wells       Wells structure.
        /// \param[in]  num_cells   Number of cells in grid.
        /// \param[in]  currentStep Number of current simulation step.
        PolymerInflowFromDeck(Opm::DeckConstPtr deck,
                              const Wells& wells,
                              const int num_cells,
                              size_t currentStep);

        /// Get inflow concentrations for all cells.
        /// \param[in]  step_start     Start of timestep.
        /// \param[in]  step_end       End of timestep.
        /// \param[out] poly_inflow_c  Injection concentrations to use for timestep, per cell.
        ///                            Must be properly sized before calling.
        virtual void getInflowValues(const double /*step_start*/,
                                     const double /*step_end*/,
                                     std::vector<double>& poly_inflow_c) const;
    private:
        SparseVector<double> sparse_inflow_;
        
        std::unordered_map<std::string, double> wellPolymerRate_;
        void setInflowValues(Opm::DeckConstPtr deck,
                             size_t currentStep);
    };


} // namespace Opm


#endif // OPM_POLYMERINFLOW_HEADER_INCLUDED
