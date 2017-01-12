/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
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

#ifndef OPM_BLACKOILDETAILS_HEADER_INCLUDED
#define OPM_BLACKOILDETAILS_HEADER_INCLUDED

#include <opm/core/linalg/ParallelIstlInformation.hpp>

namespace Opm {
namespace detail {


    inline
    std::vector<int>
    buildAllCells(const int nc)
    {
        std::vector<int> all_cells(nc);

        for (int c = 0; c < nc; ++c) { all_cells[c] = c; }

        return all_cells;
    }



    template <class PU>
    std::vector<bool>
    activePhases(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        std::vector<bool> active(maxnp, false);

        for (int p = 0; p < pu.MaxNumPhases; ++p) {
            active[ p ] = pu.phase_used[ p ] != 0;
        }

        return active;
    }



    template <class PU>
    std::vector<int>
    active2Canonical(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        std::vector<int> act2can(maxnp, -1);

        for (int phase = 0; phase < maxnp; ++phase) {
            if (pu.phase_used[ phase ]) {
                act2can[ pu.phase_pos[ phase ] ] = phase;
            }
        }

        return act2can;
    }



    inline
    double getGravity(const double* g, const int dim) {
        double grav = 0.0;
        if (g) {
            // Guard against gravity in anything but last dimension.
            for (int dd = 0; dd < dim - 1; ++dd) {
                assert(g[dd] == 0.0);
            }
            grav = g[dim - 1];
        }
        return grav;
    }


        /// \brief Compute the Euclidian norm of a vector
        /// \param it              begin iterator for the given vector
        /// \param end             end iterator for the given vector
        /// \param pinfo           In a parallel this holds the information about the data distribution.
        template <class Iterator, class Communicator>
        inline
        double euclidianNormSquared(Iterator it,
                                    const Iterator end,
                                    const Communicator& comm)
        {
            double product = 0.0 ;
            for( ; it != end; ++it ) {
                product += ( *it * *it );
            }
            return comm.sum(product);
        }

        template <class Scalar>
        inline
        double
        convergenceReduction(const std::vector< std::vector< Scalar > >& B,
                             const std::vector< std::vector< Scalar > >& tempV,
                             const std::vector< std::vector< Scalar > >& R,
                             std::vector< Scalar >& R_sum,
                             std::vector< Scalar >& maxCoeff,
                             std::vector< Scalar >& B_avg,
                             std::vector< Scalar >& maxNormWell,
                             const int nc,
                             const int np,
                             const std::vector< Scalar >& pv,
                             const std::vector< Scalar >& residual_well)
        {
            const int nw = residual_well.size() / np;
            assert(nw * np == int(residual_well.size()));

            // Do the global reductions
            {
                B_avg.resize(np);
                maxCoeff.resize(np);
                R_sum.resize(np);
                maxNormWell.resize(np);
                for ( int idx = 0; idx < np; ++idx )
                {
                    B_avg[idx] = std::accumulate( B[ idx ].begin(), B[ idx ].end(), 0.0 ) / nc;
                    R_sum[idx] = std::accumulate( R[ idx ].begin(), R[ idx ].end(), 0.0 );
                    maxCoeff[idx] = *(std::max_element( tempV[ idx ].begin(), tempV[ idx ].end() ));

                    assert(np >= np);
                    if (idx < np) {
                        maxNormWell[idx] = 0.0;
                        for ( int w = 0; w < nw; ++w ) {
                            maxNormWell[idx] = std::max(maxNormWell[idx], std::abs(residual_well[nw*idx + w]));
                        }
                    }
                }
                // Compute total pore volume
                return std::accumulate(pv.begin(), pv.end(), 0.0);
            }
        }
    } // namespace detail
} // namespace Opm

#endif // OPM_BLACKOILDETAILS_HEADER_INCLUDED
