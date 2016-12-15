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
