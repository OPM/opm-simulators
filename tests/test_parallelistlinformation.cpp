/*
  Copyright 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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

#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE OPM-ParallelIstlInformation
#include <boost/test/unit_test.hpp>
#include "DuneIstlTestHelpers.hpp"
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <functional>
#ifdef HAVE_DUNE_ISTL


template<typename T>
void runSumMaxMinTest(const T offset)
{
    const int N=100;
    int start, end, istart, iend;
    std::tie(start,istart,iend,end) = computeRegions(N);
    Opm::ParallelISTLInformation comm(MPI_COMM_WORLD);
    auto mat = create1DLaplacian(*comm.indexSet(), N, start, end, istart, iend);
    std::vector<T> x(end-start);
    assert(comm.indexSet()->size()==x.size());
    for(auto it=comm.indexSet()->begin(), itend=comm.indexSet()->end(); it!=itend; ++it)
        x[it->local()]=it->global()+offset;
    auto containers = std::make_tuple(x, x, x, x, x);
    auto operators  = std::make_tuple(Opm::Reduction::makeGlobalSumFunctor<T>(),
                                      Opm::Reduction::makeGlobalMaxFunctor<T>(),
                                      Opm::Reduction::makeGlobalMinFunctor<T>(),
                                      Opm::Reduction::makeInnerProductFunctor<T>(),
                                      Opm::Reduction::makeLInfinityNormFunctor<T>());
    auto values     = std::tuple<T,T,T,T,T>(0,0,100000, 0, 0);
    auto oldvalues  = values;
    start = offset;
    end   = start+N;
    comm.computeReduction(containers,operators,values);
    BOOST_CHECK(std::get<0>(values)==std::get<0>(oldvalues)+((N-1+2*offset)*N)/2);
    BOOST_CHECK(std::get<1>(values)==std::max(N+offset-1, std::get<1>(oldvalues)));
    BOOST_CHECK(std::get<2>(values)==std::min(offset, std::get<2>(oldvalues)));
    BOOST_CHECK(std::get<3>(values)==((end-1)*end*(2*end-1)-(start-1)*start*(2*start-1))/6+std::get<3>(oldvalues));
    // Must avoid std::abs() directly to prevent ambiguity with unsigned integers.
    Opm::Reduction::detail::MaxAbsFunctor<T> maxabsfunc;
    BOOST_CHECK(std::get<4>(values)==maxabsfunc(offset, N+offset-1));
}

BOOST_AUTO_TEST_CASE(tupleReductionTestInt)
{
    runSumMaxMinTest<int>(-200);
    runSumMaxMinTest<int>(0);
    runSumMaxMinTest<int>(20);
    runSumMaxMinTest<int>(-20);
}

BOOST_AUTO_TEST_CASE(tupleReductionTestUnsignedInt)
{
    runSumMaxMinTest<std::size_t>(0);
    runSumMaxMinTest<std::size_t>(20);
}
BOOST_AUTO_TEST_CASE(tupleReductionTestFloat)
{
    runSumMaxMinTest<float>(-200);
    runSumMaxMinTest<float>(0);
    runSumMaxMinTest<float>(20);
    runSumMaxMinTest<float>(-20);
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
    for(auto it=comm.indexSet()->begin(), itend=comm.indexSet()->end(); it!=itend; ++it)
        x[it->local()]=it->global();
    int value = 1;
    int oldvalue = value;
    comm.computeReduction(x,Opm::Reduction::makeGlobalSumFunctor<int>(),value);
    BOOST_CHECK(value==oldvalue+((N-1)*N)/2);
}
#endif
