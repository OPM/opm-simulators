/*
  Copyright 2015, 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil AS.
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


//#include <opm/autodiff/multiPhaseUpwind.hpp>
#include <algorithm>
#include <utility>


namespace Opm
{
    
    template<class Scalar, class VectorType,class IntVectorType>
    void connectionMultiPhaseUpwind(IntVectorType& upwind,
                                          const VectorType& head_diff,
                                          const VectorType& mob1,
                                          const VectorType& mob2,
                                          const Scalar transmissibility,
                                          const Scalar flux)
    {
        // Based on the paper "Upstream Differencing for Multiphase Flow in Reservoir Simulation",
        // by Yann Brenier and Jérôme Jaffré,
        // SIAM J. Numer. Anal., 28(3), 685–696.
        // DOI:10.1137/0728036
        //
        // Notation is based on this paper, except q -> flux, t -> transmissibility.

        enum { NumPhases = 3 }; // TODO: remove this restriction.

        // Get and sort the g values (also called "weights" in the paper) for this connection.
        using ValueAndIndex = std::pair<Scalar, int>;
        std::array<ValueAndIndex, NumPhases> g;
        for (int phase_idx = 0; phase_idx < NumPhases; ++phase_idx) {
            g[phase_idx] = ValueAndIndex(head_diff[phase_idx], phase_idx);
        }
        std::sort(g.begin(), g.end());

        // Compute theta and r.
        // Paper notation: subscript l -> ell (for read/searchability)
        // Note that since we index phases from 0, r is one less than in the paper.
        std::array<Scalar, NumPhases> theta;
        int r = -1;
        for (int ell = 0; ell < NumPhases; ++ell) {
            theta[ell] = flux;
            for (int j = 0; j < NumPhases; ++j) {
                if (j < ell) {
                    theta[ell] += transmissibility * (g[ell].first - g[j].first) * mob2[g[j].second];
                }
                if (j > ell) {
                    theta[ell] += transmissibility * (g[ell].first - g[j].first) * mob1[g[j].second];
                }
            }
            if (theta[ell] <= 0.0) {
                r = ell;
            } else {
                break; // r is correct, no need to continue
            }
        }

        // Set upwind array and return.
        //std::array<Scalar, NumPhases> upwind;
        for (int ell = 0; ell < NumPhases; ++ell) {
            const int phase_idx = g[ell].second;
            upwind[phase_idx] = ell > r ? 1.0 : -1.0;
        }
        //return upwind;

    }


} // namespace Opm
