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

#include <Eigen/Eigen>
#include <Eigen/Sparse>

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

        /// \brief Compute the Euclidian norm of a vector
        /// \warning In the case that num_components is greater than 1
        ///          an interleaved ordering is assumed. E.g. for each cell
        ///          all phases of that cell are stored consecutively. First
        ///          the ones for cell 0, then the ones for cell 1, ... .
        /// \param it              begin iterator for the given vector
        /// \param end             end iterator for the given vector
        /// \param num_components  number of components (i.e. phases) in the vector
        /// \param pinfo           In a parallel this holds the information about the data distribution.
        template <class Iterator>
        inline
        double euclidianNormSquared( Iterator it, const Iterator end, int num_components, const boost::any& pinfo = boost::any() )
        {
            static_cast<void>(num_components); // Suppress warning in the serial case.
            static_cast<void>(pinfo); // Suppress warning in non-MPI case.
#if HAVE_MPI
            if ( pinfo.type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>(pinfo);
                typedef typename Iterator::value_type Scalar;
                Scalar product = 0.0;
                int size_per_component = (end - it);
                size_per_component /= num_components; // two lines to supresse unused warning.
                assert((end - it) == num_components * size_per_component);

                if( num_components == 1 )
                {
                    auto component_container =
                        boost::make_iterator_range(it, end);
                    info.computeReduction(component_container,
                                           Opm::Reduction::makeInnerProductFunctor<double>(),
                                           product);
                }
                else
                {
                    auto& maskContainer = info.getOwnerMask();
                    auto mask = maskContainer.begin();
                    assert(static_cast<int>(maskContainer.size()) == size_per_component);

                    for(int cell = 0; cell < size_per_component; ++cell, ++mask)
                    {
                        Scalar cell_product = (*it) * (*it);
                        ++it;
                        for(int component=1; component < num_components;
                            ++component, ++it)
                        {
                            cell_product += (*it) * (*it);
                        }
                        product += cell_product * (*mask);
                    }
                }
                return info.communicator().sum(product);
            }
            else
#endif
            {
                double product = 0.0 ;
                for( ; it != end; ++it ) {
                    product += ( *it * *it );
                }
                return product;
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
