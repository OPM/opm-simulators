/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELLDENSITYSEGMENTED_HEADER_INCLUDED
#define OPM_WELLDENSITYSEGMENTED_HEADER_INCLUDED

#include <vector>

struct Wells;

namespace Opm
{

    class WellStateFullyImplicitBlackoil;
    class WellStateFullyImplicitBlackoilSolvent;
    struct PhaseUsage;


    /// A class giving a well model, by which we mean a way to compute
    /// the pressure deltas of each perforation and the bottom-hole
    /// pressure. This class contains an explicit model, that uses a
    /// different density for each well segment, that is between each
    /// pair of perforations.
    class WellDensitySegmented
    {
    public:
        /// Compute well segment densities
        /// Notation: N = number of perforations, P = number of phases.
        /// \param[in] wells        struct with static well info
        /// \param[in] wstate       dynamic well solution information, only perfRates() is used
        /// \param[in] phase_usage  specifies which phases are active and not
        /// \param[in] b_perf       inverse ('little b') formation volume factor, size NP, P values per perforation
        /// \param[in] rsmax_perf   saturation point for rs (gas in oil) at each perforation, size N
        /// \param[in] rvmax_perf   saturation point for rv (oil in gas) at each perforation, size N
        /// \param[in] surf_dens    surface densities for active components, size NP, P values per perforation
        static std::vector<double> computeConnectionDensities(const Wells& wells,
                                                              const WellStateFullyImplicitBlackoil& wstate,
                                                              const PhaseUsage& phase_usage,
                                                              const std::vector<double>& b_perf,
                                                              const std::vector<double>& rsmax_perf,
                                                              const std::vector<double>& rvmax_perf,
                                                              const std::vector<double>& surf_dens_perf);



        /// Compute well segment densities for solvent model
        /// Notation: N = number of perforations, P = number of phases.
        /// \param[in] wells        struct with static well info
        /// \param[in] wstate       dynamic well solution information, perfRates() and solventFraction() is used
        /// \param[in] phase_usage  specifies which phases are active and not
        /// \param[in] b_perf       inverse ('little b') formation volume factor, size NP, P values per perforation
        /// \param[in] rsmax_perf   saturation point for rs (gas in oil) at each perforation, size N
        /// \param[in] rvmax_perf   saturation point for rv (oil in gas) at each perforation, size N
        /// \param[in] surf_dens    surface densities for active components, size NP, P values per perforation
        static std::vector<double> computeConnectionDensities(const Wells& wells,
                                                              const WellStateFullyImplicitBlackoilSolvent& wstate,
                                                              const PhaseUsage& phase_usage,
                                                              const std::vector<double>& b_perf,
                                                              const std::vector<double>& rsmax_perf,
                                                              const std::vector<double>& rvmax_perf,
                                                              const std::vector<double>& surf_dens_perf);




        /// Compute pressure deltas.
        /// Notation: N = number of perforations, P = number of phases.
        /// \param[in] wells        struct with static well info
        /// \param[in] z_perf       depth values for each perforation, size N
        /// \param[in] dens_perf    densities for each perforation, size N (typically computed using computeConnectionDensities)
        /// \param[in] gravity      gravity acceleration constant
        static std::vector<double> computeConnectionPressureDelta(const Wells& wells,
                                                                  const std::vector<double>& z_perf,
                                                                  const std::vector<double>& dens_perf,
                                                                  const double gravity);
    };

} // namespace Opm

#endif // OPM_WELLDENSITYSEGMENTED_HEADER_INCLUDED
