/*
  Copyright 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU

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

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE OPM-ParallelIstlInformation
#include <boost/test/unit_test.hpp>
#if HAVE_MPI
#include <mpi.h>
#else
#error "This file needs to compiled with MPI support!"
#endif
#include "DuneIstlTestHelpers.hpp"
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <functional>
#ifdef HAVE_DUNE_ISTL
BOOST_AUTO_TEST_CASE(tupleReductionTest)
{
    int N=100;
    int start, end, istart, iend;
    std::tie(start,istart,iend,end) = computeRegions(N);
    Opm::ParallelISTLInformation comm(MPI_COMM_WORLD);
    auto mat = create1DLaplacian(*comm.indexSet(), N, start, end, istart, iend);
    std::vector<int> x(end-start);
    assert(comm.indexSet()->size()==x.size());
    for(auto i=comm.indexSet()->begin(), iend=comm.indexSet()->end(); i!=iend; ++i)
        x[i->local()]=i->global();
    auto containers = std::make_tuple(x, x, x);
    auto operators  = std::make_tuple(Opm::Reduction::MaskIDOperator<std::plus<int> >(),
                                      Opm::Reduction::MaskToMinOperator<std::greater<int> >(),
                                      Opm::Reduction::MaskToMaxOperator<std::less< int>  >());
    auto values     = std::make_tuple(0,0,100000);
    auto oldvalues  = values;
    comm.computeReduction(containers,operators,values);
    BOOST_CHECK(std::get<0>(values)==std::get<0>(oldvalues)+((N-1)*N)/2);
    BOOST_CHECK(std::get<1>(values)==std::min(0, std::get<1>(oldvalues)));
    BOOST_CHECK(std::get<2>(values)==std::max(N, std::get<2>(oldvalues)));
}
BOOST_AUTO_TEST_CASE(singleContainerReductionTest)
{
    int N=100;
    int start, end, istart, iend;
    std::tie(start,istart,iend,end) = computeRegions(N);
    Opm::ParallelISTLInformation comm(MPI_COMM_WORLD);
    auto mat = create1DLaplacian(*comm.indexSet(), N, start, end, istart, iend);
    std::vector<int> x(end-start);
    assert(comm.indexSet()->size()==x.size());
    for(auto i=comm.indexSet()->begin(), iend=comm.indexSet()->end(); i!=iend; ++i)
        x[i->local()]=i->global();
    auto containers = std::make_tuple(x, x, x);
    auto operators  = std::make_tuple(Opm::Reduction::MaskIDOperator<std::plus<int> >(),
                                      Opm::Reduction::MaskToMinOperator<std::greater<int> >(),
                                      Opm::Reduction::MaskToMaxOperator<std::less< int>  >());
    int value = 1;
    int oldvalue = value;
    comm.computeReduction(x,Opm::Reduction::MaskIDOperator<std::plus<int> >(),value);
    BOOST_CHECK(value==oldvalue+((N-1)*N)/2);
}
#endif
