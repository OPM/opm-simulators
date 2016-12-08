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

#ifndef OPM_BLACKOILLEGACYDETAILS_HEADER_INCLUDED
#define OPM_BLACKOILLEGACYDETAILS_HEADER_INCLUDED

#include <opm/core/linalg/ParallelIstlInformation.hpp>

namespace Opm {
namespace detail {
        /// \brief Compute the L-infinity norm of a vector
        /// \warn This function is not suitable to compute on the well equations.
        /// \param a The container to compute the infinity norm on.
        ///          It has to have one entry for each cell.
        /// \param info In a parallel this holds the information about the data distribution.
        template <class ADB>
        inline
        double infinityNorm( const ADB& a, const boost::any& pinfo = boost::any() )
        {
            static_cast<void>(pinfo); // Suppress warning in non-MPI case.
#if HAVE_MPI
            if ( pinfo.type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& real_info =
                    boost::any_cast<const ParallelISTLInformation&>(pinfo);
                double result=0;
                real_info.computeReduction(a.value(), Reduction::makeLInfinityNormFunctor<double>(), result);
                return result;
            }
            else
#endif
            {
                if( a.value().size() > 0 ) {
                    return a.value().matrix().template lpNorm<Eigen::Infinity> ();
                }
                else { // this situation can occur when no wells are present
                    return 0.0;
                }
            }
        }

        /// \brief Compute the reduction within the convergence check.
        /// \param[in] B     A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. B.col(i) contains the values
        ///                  for phase i.
        /// \param[in] tempV A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. tempV.col(i) contains the
        ///                   values
        ///                  for phase i.
        /// \param[in] R     A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. B.col(i) contains the values
        ///                  for phase i.
        /// \param[out] R_sum An array of size MaxNumPhases where entry i contains the sum
        ///                   of R for the phase i.
        /// \param[out] maxCoeff An array of size MaxNumPhases where entry i contains the
        ///                   maximum of tempV for the phase i.
        /// \param[out] B_avg An array of size MaxNumPhases where entry i contains the average
        ///                   of B for the phase i.
        /// \param[out] maxNormWell The maximum of the well flux equations for each phase.
        /// \param[in]  nc    The number of cells of the local grid.
        /// \return The total pore volume over all cells.
        inline
        double
        convergenceReduction(const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& B,
                             const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& tempV,
                             const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& R,
                             std::vector<double>& R_sum,
                             std::vector<double>& maxCoeff,
                             std::vector<double>& B_avg,
                             std::vector<double>& maxNormWell,
                             int nc,
                             int np,
                             const std::vector<double> pv,
                             std::vector<double> residual_well)
        {
            const int nw = residual_well.size() / np;
            assert(nw * np == int(residual_well.size()));

            // Do the global reductions
#if 0 // HAVE_MPI
            if ( linsolver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>(linsolver_.parallelInformation());

                // Compute the global number of cells and porevolume
                std::vector<int> v(nc, 1);
                auto nc_and_pv = std::tuple<int, double>(0, 0.0);
                auto nc_and_pv_operators = std::make_tuple(Opm::Reduction::makeGlobalSumFunctor<int>(),
                                                           Opm::Reduction::makeGlobalSumFunctor<double>());
                auto nc_and_pv_containers  = std::make_tuple(v, pv);
                info.computeReduction(nc_and_pv_containers, nc_and_pv_operators, nc_and_pv);

                for ( int idx = 0; idx < np; ++idx )
                {
                    auto values     = std::tuple<double,double,double>(0.0 ,0.0 ,0.0);
                    auto containers = std::make_tuple(B.col(idx),
                                                      tempV.col(idx),
                                                      R.col(idx));
                    auto operators  = std::make_tuple(Opm::Reduction::makeGlobalSumFunctor<double>(),
                                                      Opm::Reduction::makeGlobalMaxFunctor<double>(),
                                                      Opm::Reduction::makeGlobalSumFunctor<double>());
                    info.computeReduction(containers, operators, values);
                    B_avg[idx]       = std::get<0>(values)/std::get<0>(nc_and_pv);
                    maxCoeff[idx]    = std::get<1>(values);
                    R_sum[idx]       = std::get<2>(values);
                    assert(np >= np);
                    if (idx < np) {
                        maxNormWell[idx] = 0.0;
                        for ( int w = 0; w < nw; ++w ) {
                            maxNormWell[idx]  = std::max(maxNormWell[idx], std::abs(residual_well[nw*idx + w]));
                        }
                    }
                }
                info.communicator().max(maxNormWell.data(), np);
                // Compute pore volume
                return std::get<1>(nc_and_pv);
            }
            else
#endif
            {
                B_avg.resize(np);
                maxCoeff.resize(np);
                R_sum.resize(np);
                maxNormWell.resize(np);
                for ( int idx = 0; idx < np; ++idx )
                {
                    B_avg[idx] = B.col(idx).sum()/nc;
                    maxCoeff[idx] = tempV.col(idx).maxCoeff();
                    R_sum[idx] = R.col(idx).sum();

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

        /// \brief Compute the L-infinity norm of a vector representing a well equation.
        /// \param a The container to compute the infinity norm on.
        /// \param info In a parallel this holds the information about the data distribution.
        template <class ADB>
        inline
        double infinityNormWell( const ADB& a, const boost::any& pinfo )
        {
            static_cast<void>(pinfo); // Suppress warning in non-MPI case.
            double result=0;
            if( a.value().size() > 0 ) {
                result = a.value().matrix().template lpNorm<Eigen::Infinity> ();
            }
#if HAVE_MPI
            if ( pinfo.type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& real_info =
                    boost::any_cast<const ParallelISTLInformation&>(pinfo);
                result = real_info.communicator().max(result);
            }
#endif
            return result;
        }
    } // namespace detail
} // namespace Opm

#endif // OPM_BLACKOILDETAILS_HEADER_INCLUDED
